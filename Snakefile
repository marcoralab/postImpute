'''Snakefile for MIS post-imputation QC
   Version 0.1.1'''
try:
    from scripts.parse_config import parser
except:
    from workflow.scripts.parse_config import parser
from getpass import getuser

import pandas as pd
import os
import socket

RWD = os.getcwd()

configfile: "config/config.yaml"
SAMPLES = pd.DataFrame.from_records(config["SAMPLES"], index = "COHORT")

BPLINK = ["bed", "bim", "fam"]

CHROM, SAMPLE, INPATH, KEEP_COMMAND = parser(config)

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
    sampfilt += com["bcftools"] + " view --samples-file "
    if config["exclude_samp"]:
        sampfilt += "^{}".format(config["exclude_samp"])
    if config["include_samp"]:
        sampfilt += "{}".format(config["include_samp"])

plink_bycohort = "data/{sample}_chrall_filtered.{ext}"
plink_merged = "data/all_chrall_filtered.{ext}"

outs = dict(
    stat_report=expand("stats/{sample}_impStats.html", sample=SAMPLE),
    vcf_bycohort=expand("data/{sample}_chrall_filtered.vcf.gz", sample=SAMPLE),
    vcf_merged="data/all_chrall_filtered.vcf.gz",
    bgen_bycohort=expand("data/{sample}_chrall_filtered.bgen", sample=SAMPLE),
    bgen_merged="data/merged/merged_chrall_filtered.bgen",
    plink_bycohort=expand(plink_bycohort, sample=SAMPLE, ext=BPLINK),
    plink_merged=expand(plink_merged, ext=BPLINK)
            )

outputs = [outs[x] for x in config["outputs"]]
outputs = flatten(outputs)

rule all:
    input: outputs


rule unzip:
    input: INPATH + "{sample}/chr_{chrom}.zip"
    output:
        vcf = "input/{sample}/chr{chrom}.dose.vcf.gz",
        info = "input/{sample}/chr{chrom}.info.gz"
    params:
        passwd = lambda wildcards: SAMPLES.loc[wildcards.sample]['JOB']['pwd'],
        odir = "input/{sample}"
    conda: 'workflow/envs/p7z.yaml'
    shell:
        r'''
7za e {input} -p{params.passwd} -o{params.odir}
for FILE in chunks-excluded.txt snps-excluded.txt typed-only.txt chr_{{1..22}}.log; do
  INFILE="{INPATH}{wildcards.sample}/$FILE"
  OUTFILE="{params.odir}/$FILE"
  if test -f "$INFILE"; then
    echo copying $INFILE to $OUTFILE 
    cp $INFILE $OUTFILE
  fi
done
'''

rule stats:
    input:
        markdown = "workflow/scripts/Post_imputation.Rmd",
        info = expand("input/{sample}/chr{chrom}.info.gz", sample=SAMPLE, chrom=CHROM)
    output:
        outfile = "stats/{sample}_impStats.html"
    params:
        rwd = RWD,
        path = "input/{sample}/",
        chrom = config["chroms"],
        cohort = "{sample}",
        maf = config["qc"]["maf"],
        rsq = config["qc"]["rsq"],
        rsq2 = config["qc"]["rsq2"],
        sampsize = config["qc"]["sampsize"],
        out = "{sample}_impStats.html",
        output_dir = "stats"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/RenderPostImputationReport.R"

# Sample filtering rules

rule indexinitial:
    input: rules.unzip.output.vcf
    output: "input/{sample}/chr{chrom}.dose.vcf.gz.tbi"
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools index -t {input}"

rule fixheaders:
    input:
        vcf = rules.unzip.output.vcf,
        tbi = rules.indexinitial.output,
    output:
        vcf = temp("temp/fixedheader/{sample}/chr{chrom}.dose.vcf.gz"),
        tbi = temp("temp/fixedheader/{sample}/chr{chrom}.dose.vcf.gz.tbi"),
    threads: 1
    conda: "workflow/envs/bcftools.yaml"
    shell:
        r"""
if $(zcat {input} | head -n 40 | grep -q "##FILTER"); then
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

filter_annotate = ("bcftools annotate -i \"%FILTER='GENOTYPED' || "
                   "{params.filt}\" -Oz -o {output} "
                   "--set-id '%CHROM:%POS:%REF:%ALT' --threads 8")
filter_out = "data/by_chrom/{sample}_chr{chrom}_filtered.vcf.gz"

if sampfilt:
    rule filters:
        input:
            vcf = rules.fixheaders.output.vcf,
            tbi = rules.fixheaders.output.tbi
        output: temp(filter_out)
        params:
            filt = qualfilt,
            sf = sampfilt
        threads: 8
        conda: "workflow/envs/bcftools.yaml"
        shell:
            "{params.sf} --force-samples -Oz --threads 8 {input.vcf} | "
            "" + filter_annotate

else:
    rule filters:
        input:
            vcf = rules.fixheaders.output.vcf,
            tbi = rules.fixheaders.output.tbi
        output: temp(filter_out)
        params:
            filt = qualfilt,
            sf = sampfilt
        threads: 8
        conda: "workflow/envs/bcftools.yaml"
        shell: filter_annotate + " {input.vcf}"

# defaults for renaming:
renamed_cat = "data/by_chrom/{{sample}}_chr{chrom}_filtered.vcf.gz"
renamed = "data/by_chrom/{sample}_chr{chrom}_filtered.vcf.gz"
renamed_merge = "data/by_chrom/{sample}_chr{{chrom}}_filtered.vcf.gz"
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
        renamed_cat = "data/by_chrom/{{sample}}_chr{chrom}_filtered_fixedIDs.vcf.gz"
        renamed = "data/by_chrom/{sample}_chr{chrom}_filtered_fixedIDs.vcf.gz"
        renamed_merge = "data/by_chrom/{sample}_chr{{chrom}}_filtered_fixedIDs.vcf.gz"
    elif config['rename'] and not type(config['rename']) is dict:
        print("Manualy renameing samples")
        renamefile = config['rename']
        renamed_cat = "data/by_chrom/{{sample}}_chr{chrom}_filtered_renamed.vcf.gz"
        renamed = "data/by_chrom/{sample}_chr{chrom}_filtered_renamed.vcf.gz"
        renamed_merge = "data/by_chrom/{sample}_chr{{chrom}}_filtered_renamed.vcf.gz"

rule rename:
    input:
        vcf = rules.filters.output,
        mapping = renamefile if rename_tf else "/dev/null"
    output: temp("data/by_chrom/{sample}_chr{chrom}_filtered_renamed.vcf.gz")
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools reheader --samples {input.mapping} -o {output} -Oz {input.vcf}"

rule fixHeader:
    input:
        vcf = rules.filters.output,
        mapping = automap if automap_tf else "/dev/null"
    output:
        mapping = "data/by_chrom/{sample}_chr{chrom}_vcfmap.tsv",
        reheader = temp("data/by_chrom/{sample}_chr{chrom}_vcfreheader.txt")
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/fix_HRCvcf.R"

rule renameAuto:
    input:
        vcf = rules.filters.output,
        header = rules.fixHeader.output.reheader
    output:
        fixed = temp("data/by_chrom/{sample}_chr{chrom}_filtered_fixedIDs.vcf.gz"),
    conda: "workflow/envs/bcftools.yaml"
    shell:
        "bcftools reheader --samples {output.reheader} {input.vcf} | "
        "bcftools view -o {output.fixed} -Oz"

rule concat_chroms_samp:
    input: expand(renamed_cat, chrom=CHROM)
    output: "data/{sample}_chrall_filtered.vcf.gz"
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
        vcf = expand(renamed_merge, sample=SAMPLE),
        tbi = expand(renamed_merge + ".tbi", sample=SAMPLE)
    output: "data/by_chrom/all_chr{chrom}_filtered.vcf.gz"
    threads: 8
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools merge -m none -o {output} -Oz --threads 8 {input.vcf}"

rule concat_chroms_all:
    input: expand("data/by_chrom/all_chr{chrom}_filtered.vcf.gz", chrom=CHROM)
    output: "data/all_chrall_filtered.vcf.gz"
    threads: 8
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools concat -o {output} -Oz --threads 8 {input}"

rule make_plink_all:
    input: rules.concat_chroms_all.output
    output: expand("data/all_chrall_filtered.{ext}", ext=BPLINK)
    params:
        out_plink = "data/all_chrall_filtered",
        ID = "--id-delim" if automap_tf else "--double-id"
    threads: 10
    conda: "workflow/envs/plink.yaml"
    shell:
        "plink --keep-allele-order --vcf {input} {params.ID} --memory 10000 --threads 10 --make-bed "
        "--out {params.out_plink}"

rule make_plink_samp:
    input: rules.concat_chroms_samp.output
    output: expand("data/{{sample}}_chrall_filtered.{ext}", ext=BPLINK)
    params:
        out_plink = "data/{sample}_chrall_filtered",
        ID = "--id-delim" if automap_tf else "--double-id "
    threads: 10
    conda: "workflow/envs/plink.yaml"
    shell:
        "plink --keep-allele-order --vcf {input} {params.ID} --memory 10000 --threads 10 --make-bed "
        "--out {params.out_plink}"

# If bgen outputs are requested
include: "workflow/rules/bgen.smk"
