'''Snakefile for MIS post-imputation QC
   Version 0.1.1'''

from workflow.scripts.parse_config import parser_postImpute
from getpass import getuser

import pandas as pd
import os
import socket

RWD = os.getcwd()

# import ipdb; ipdb.set_trace()

configfile: "config/config.yaml"
if 'SAMPLES' in config:
    zipped = True
    S2 = [{'COHORT': x, 'JOB': y} for x, y in config["SAMPLES"].items()]
    SAMPLES = pd.DataFrame.from_records(S2, index = "COHORT")
else:
    zipped = False

BPLINK = ["bed", "bim", "fam"]

CHROM, COHORT, INPATH, KEEP_COMMAND = parser_postImpute(config)
fixheaders = config['fixheaders'] if 'fixheaders' in config else True

#import ipdb; ipdb.set_trace()

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

# remove subjects if sample filtering file is provided
sampfilt = ""
if config["exclude_samp"] or config["include_samp"]:
    sampfilt += "bcftools view --samples-file "
    if config["exclude_samp"]:
        sampfilt += "^{}".format(config["exclude_samp"])
    if config["include_samp"]:
        sampfilt += "{}".format(config["include_samp"])

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
# import ipdb; ipdb.set_trace()
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
    return expand("{impute_dir}/input/{cohort}/chr{chrom}.info.gz",
        impute_dir=wildcards["impute_dir"],
        cohort=wildcards["cohort"],
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
