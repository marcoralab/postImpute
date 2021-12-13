'''Snakefile for MIS post-imputation QC'''

from workflow.scripts.parse_config import parser_postImpute
from getpass import getuser

import pandas as pd
import os
import socket
import re
import glob
import shutil


RWD = os.getcwd()

configfile: "config/config.yaml"

zipped = 'SAMPLES' in config else

BPLINK = ["bed", "bim", "fam"]

# --- Process chromosomes config ---


def parse_chrom(chrs):
    clist = [x.split(":") for x in chrs.split(",")]
    parsed = []
    for chrs in clist:
        if len(chrs) == 2:
            chrs = [str(c) for c in range(int(chrs[0]), int(chrs[1]) + 1)]
        elif len(chrs) != 1:
            raise ValueError("Invalid chromosome list.")
        parsed += chrs
    return parsed


if 'chrom' not in config:
    CHROM = parse_chrom('1:22')
else:
    CHROM = parse_chrom(config['chroms'])

# --- Process keep and remove ---


def build_sampfilt_vcf(config):
    def filtstr(x):
        if not (x in config and config[x] and config[x] is not None):
            return None
        elif os.path.isfile(config[x]):
            paramstr = '' if x == 'include_samp' else '^'
            return '{}{}'.format(paramstr, os.path.normpath(config[x]))
        else:
            filt = 'inclusion' if x == 'include_samp' else 'exclusion'
            raise Exception("Invalid {} list: {}.".format(filt, config[x]))
    x = [filtstr(x) for x in ['include_samp', 'exclude_samp']]
    x = [i for i in x if i is not None]
    if len(x) > 1:
        raise Exception('Cannot have both sample removal and exclusion.')
    return 'bcftools view --samples-file '+ x[0] if x else ''

sampfilt = build_sampfilt_vcf(config)

# --- Process input files from config ---


def build_samp(in_path, samples = None):
    p_abs = os.path.abspath(in_path)
    p = [directory for directory,y,files in os.walk(p_abs)
         if any(".zip" in f for f in files)]
    if len(p) == 0:
        p = [directory for directory,y,files in os.walk(p_abs)
             if any("info.gz" in f for f in files)]
    if samples is not None:
        p = [x for x in p if x in samples]
    return [os.path.basename(x) for x in p]


INPATH = os.path.abspath(config["directory"]) + '/' # normalize input dir

if "SAMPLES" in config:
    COHORT = build_samp(INPATH, [*config["SAMPLES"]])
else:
    COHORT = build_samp(INPATH)

# --- Done processing ---

fixheaders = config['fixheaders'] if 'fixheaders' in config else True


def flatten(nested):
    flat = []
    for el in nested:
        if not isinstance(el, list):
            flat.append(el)
        else:
            flat += flatten(el)
    return flat


# Change filtering based on presence of rsq2
qualfilt = "(R2 >= {R2} && MAF >= {MAF})".format(
    R2=config["qc"]["rsq"], MAF=config["qc"]["maf"])
if config["qc"]["rsq2"] and config["qc"]["rsq2"] != 'NA':
    qualfilt += " || (R2 >= {R2} && MAF < {MAF})".format(
        R2=config["qc"]["rsq2"], MAF=config["qc"]["maf"])

plink_bycohort = "{impute_dir}/data/{cohort}_chrall_filtered.{ext}"
plink_merged = "{impute_dir}/data/all_chrall_filtered.{ext}"

outs = dict(
    stat_report=expand("{impute_dir}/stats/{cohort}_impStats.html", cohort=COHORT, impute_dir=config['directory']),
    vcf_bycohort=expand("{impute_dir}/data/{cohort}_chrall_filtered.vcf.gz", cohort=COHORT, impute_dir=config['directory']),
    vcf_merged=expand("{impute_dir}/data/all_chrall_filtered.vcf.gz", impute_dir=config['directory']),
    bgen_bycohort=expand("{impute_dir}/data/{cohort}_chrall_filtered.bgen", cohort=COHORT, impute_dir=config['directory']),
    bgen_merged=expand("{impute_dir}/data/merged/merged_chrall_filtered.bgen", impute_dir=config['directory']),
    plink_bycohort=expand(plink_bycohort, cohort=COHORT, ext=BPLINK, impute_dir=config['directory']),
    plink_merged=expand(plink_merged, ext=BPLINK, impute_dir=config['directory'])
            )

outputs = [outs[x] for x in config["outputs"]]
outputs = flatten(outputs)

rule all:
    input: outputs

if zipped:
    rule unzip:
        input: "{cohort}/chr_{chrom}.zip"
        output:
            vcf = "{impute_dir}/input/{cohort}/chr{chrom}.dose.vcf.gz",
            info = "{impute_dir}/input/{cohort}/chr{chrom}.info.gz"
        params:
            passwd = lambda wildcards: SAMPLES.loc[wildcards.cohort]['JOB']['pwd'],
            odir = "{impute_dir}/input/{cohort}",
            INPATH = INPATH
        conda: 'workflow/envs/p7z.yaml'
        shell:
            r'''
            7za e {input} -p{params.passwd} -o{params.odir}
            for FILE in chunks-excluded.txt snps-excluded.txt typed-only.txt chr_{{1..22}}.log; do
              INFILE="{params.INPATH}/{wildcards.cohort}/$FILE"
              OUTFILE="{params.odir}/$FILE"
              if test -f "$INFILE"; then
                echo copying $INFILE to $OUTFILE
                cp $INFILE $OUTFILE
              fi
            done
            '''


def stats_input(wildcards):
    if zipped:
        return expand("{impute_dir}/input/{{cohort}}/chr{chrom}.info.gz",
            impute_dir=wildcards["impute_dir"],
            chrom=CHROM)
    else:
        return expand(INPATH + "{{cohort}}/chr{chrom}.info.gz",
            chrom=CHROM)


rule stats:
    input:
        info = stats_input
    output:
        outfile = "{impute_dir}/stats/{cohort}_impStats.html"
    params:
        rwd = RWD,
        path = "{impute_dir}/input/{cohort}/",
        chrom = CHROM,
        cohort = "{cohort}",
        maf = config["qc"]["maf"],
        rsq = config["qc"]["rsq"],
        rsq2 = config["qc"]["rsq2"],
        sampsize = config["qc"]["sampsize"],
        out = "{cohort}_impStats.html",
        output_dir = "{impute_dir}/stats",
        reading = "P"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/Post_imputation.Rmd"

# Sample filtering rules
startfile =  "{impute_dir}/input/{cohort}/chr{chrom}.dose.vcf.gz" if zipped else INPATH + "{cohort}/chr{chrom}.dose.vcf.gz"

rule indexinitial:
    input: startfile
    output: startfile + ".tbi"
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools index -t {input}"

rule fixheaders:
    input:
        vcf = startfile,
        tbi = rules.indexinitial.output,
    output:
        vcf = temp("{impute_dir}/temp/fixedheader/{cohort}/chr{chrom}.dose.vcf.gz"),
        tbi = temp("{impute_dir}/temp/fixedheader/{cohort}/chr{chrom}.dose.vcf.gz.tbi"),
    threads: 1
    conda: "workflow/envs/bcftools.yaml"
    shell:
        r"""
if $(zcat {input} | head -n 4000 | grep -q -e "##FILTER" -e "##INFO=<ID=IMPUTED"); then
  cp {input.vcf} {output.vcf}
  cp {input.tbi} {output.tbi}
else
  zcat {input.vcf} | sed '/contig/ a\
##FILTER=<ID=GENOTYPED,Description="Marker was genotyped AND imputed">\
##FILTER=<ID=GENOTYPED_ONLY,Description="Marker was genotyped but NOT imputed">' | \
bcftools view -Oz -o {output.vcf}
  bcftools index -t {output.vcf}
fi
"""

filter_out = "{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered.vcf.gz"


if sampfilt:
    rule filters:
        input:
            vcf = rules.fixheaders.output.vcf if fixheaders else startfile,
            tbi = rules.fixheaders.output.tbi if fixheaders else startfile + '.tbi'
        output: temp(filter_out)
        params:
            filt = qualfilt,
            sf = sampfilt
        threads: 8
        conda: "workflow/envs/bcftools.yaml"
        shell:
            r'''
if zcat {input.vcf} | head -n 4000 | grep -q "##INFO=<ID=IMPUTED"; then
  {params.sf} --force-samples -Oz --threads 8 {input.vcf} | \
  bcftools annotate -i "INFO/TYPED=1 || INFO/TYPED_ONLY=1 || {params.filt}" -Oz -o {output} --set-id '%CHROM:%POS:%REF:%ALT' --threads 8
else
  {params.sf} --force-samples -Oz --threads 8 {input.vcf} | \
  bcftools annotate -i "%FILTER='GENOTYPED' || %FILTER='GENOTYPED_ONLY' || {params.filt}" -Oz -o {output} --set-id '%CHROM:%POS:%REF:%ALT' --threads 8
fi'''

else:
    rule filters:
        input:
            vcf = rules.fixheaders.output.vcf if fixheaders else startfile,
            tbi = rules.fixheaders.output.tbi if fixheaders else startfile + '.tbi'
        output: temp(filter_out)
        params:
            filt = qualfilt
        threads: 8
        conda: "workflow/envs/bcftools.yaml"
        shell:
            r'''
if zcat {input.vcf} | head -n 4000 | grep -q "##INFO=<ID=IMPUTED"; then
  echo "Processing Minimac4 file"
  bcftools annotate -i "INFO/TYPED=1 || INFO/TYPED_ONLY=1 || {params.filt}" -Oz -o {output} --set-id '%CHROM:%POS:%REF:%ALT' --threads 8 {input.vcf}
else
  echo "Processing Minimac3 file"
  bcftools annotate -i "%FILTER='GENOTYPED' || %FILTER='GENOTYPED_ONLY' || {params.filt}" -Oz -o {output} --set-id '%CHROM:%POS:%REF:%ALT' --threads 8 {input.vcf}
fi'''

# defaults for renaming:
renamed_cat = "{{impute_dir}}/data/by_chrom/{{cohort}}_chr{chrom}_filtered.vcf.gz"
renamed = "{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered.vcf.gz"
renamed_merge = "{impute_dir}/data/by_chrom/{cohort}_chr{{chrom}}_filtered.vcf.gz"
automap_tf = False
rename_tf = False

#override defaults:
if 'rename' in config:
    if (config['rename'] is not None and
        'automap' in config['rename'] and
        config['rename']['automap']):
        print("Automatically renaming samples.")
        automap_tf = True
        automap = config['rename']['automap']
        renamed_cat = "{{impute_dir}}/data/by_chrom/{{cohort}}_chr{chrom}_filtered_fixedIDs.vcf.gz"
        renamed = "{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered_fixedIDs.vcf.gz"
        renamed_merge = "{impute_dir}/data/by_chrom/{cohort}_chr{{chrom}}_filtered_fixedIDs.vcf.gz"
    elif config['rename'] and not type(config['rename']) is dict:
        print("Manualy renameing samples")
        renamefile = config['rename']
        renamed_cat = "{{impute_dir}}/data/by_chrom/{{cohort}}_chr{chrom}_filtered_renamed.vcf.gz"
        renamed = "{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered_renamed.vcf.gz"
        renamed_merge = "{impute_dir}/data/by_chrom/{cohort}_chr{{chrom}}_filtered_renamed.vcf.gz"

rule rename:
    input:
        vcf = rules.filters.output,
        mapping = renamefile if rename_tf else "/dev/null"
    output: temp("{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered_renamed.vcf.gz")
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools reheader --samples {input.mapping} -o {output} -Oz {input.vcf}"

rule fixHeader:
    input:
        vcf = rules.filters.output,
        mapping = automap if automap_tf else "/dev/null"
    output:
        mapping = "{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_vcfmap.tsv",
        reheader = temp("{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_vcfreheader.txt")
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/fix_HRCvcf.R"

rule renameAuto:
    input:
        vcf = rules.filters.output,
        header = rules.fixHeader.output.reheader
    output:
        fixed = temp("{impute_dir}/data/by_chrom/{cohort}_chr{chrom}_filtered_fixedIDs.vcf.gz"),
    conda: "workflow/envs/bcftools.yaml"
    shell:
        "bcftools reheader --samples {output.reheader} {input.vcf} | "
        "bcftools view -o {output.fixed} -Oz"

rule concat_chroms_samp:
    input: expand(renamed_cat, chrom=CHROM)
    output: "{impute_dir}/data/{cohort}_chrall_filtered.vcf.gz"
    threads: 8
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools concat -o {output} -Oz --threads 8 {input}"

rule index_samples_chrom:
    input: renamed
    output: renamed + ".tbi"
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools index -t {input}"

rule merge_samples_chrom:
    input:
        vcf = expand(renamed_merge, cohort=COHORT, impute_dir=config['directory']),
        tbi = expand(renamed_merge + ".tbi", cohort=COHORT, impute_dir=config['directory'])
    output: "{impute_dir}/data/by_chrom/all_chr{chrom}_filtered.vcf.gz"
    threads: 8
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools merge -m none -o {output} -Oz --threads 8 {input.vcf}"

rule concat_chroms_all:
    input: expand("data/by_chrom/all_chr{chrom}_filtered.vcf.gz", chrom=CHROM, impute_dir=config['directory'])
    output: "{impute_dir}/data/all_chrall_filtered.vcf.gz"
    threads: 8
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools concat -o {output} -Oz --threads 8 {input}"

rule make_plink_all:
    input: rules.concat_chroms_all.output
    output: multiext("{impute_dir}/data/all_chrall_filtered", ".bed", ".bim", ".fam")
    params:
        out_plink = "{impute_dir}/data/all_chrall_filtered",
        ID = "--id-delim" if automap_tf else "--double-id"
    threads: 10
    conda: "workflow/envs/plink.yaml"
    shell:
        "plink --keep-allele-order --vcf {input} {params.ID} --memory 10000 --threads 10 --make-bed "
        "--out {params.out_plink}"

rule make_plink_samp:
    input: rules.concat_chroms_samp.output
    output: multiext("{impute_dir}/data/{cohort}_chrall_filtered", ".bed", ".bim", ".fam")
    params:
        out_plink = "{impute_dir}/data/{cohort}_chrall_filtered",
        ID = "--id-delim" if automap_tf else "--double-id "
    threads: 10
    conda: "workflow/envs/plink.yaml"
    shell:
        "plink --keep-allele-order --vcf {input} {params.ID} --memory 10000 --threads 10 --make-bed "
        "--out {params.out_plink}"

# If bgen outputs are requested
include: "workflow/rules/bgen.smk"
