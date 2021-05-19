# ```snakemake file to filter vcf files```

import pandas as pd

configfile: "config/config_posthrc.yaml"
SAMPLES = pd.DataFrame.from_records(config["SAMPLES"], index = "COHORT")
CHR = range(1,23)

rule all:
    input:
        expand("input/{sample}/chr{chr}.dose.vcf.gz", sample = SAMPLES.index.tolist(), chr = CHR),
        expand("input/{sample}/chr{chr}.info.gz", sample = SAMPLES.index.tolist(), chr = CHR),
        # expand("temp/{sample}/chr{chr}.dose.vcf.gz", sample = SAMPLE, chr = CHR),
        # expand("temp/{sample}/chrall.dose.vcf.gz", sample = SAMPLE),
        # expand("HRC_impute/final/{sample}_chrall{bplink}", sample = SAMPLES.index.tolist(), bplink = ['.bed', '.bim', '.fam']),
        # expand("HRC_impute/final/{sample}_chrall_fixed.fam", sample = SAMPLES.index.tolist(), bplink = ['.bed', '.bim', '.fam']),
        # expand("HRC_impute/final/{sample}_snp_info.txt", sample = SAMPLES.index.tolist())


rule unzip:
    input: "HRC_impute/{sample}/chr_{chr}.zip"
    output:
        vcf = "input/{sample}/chr{chr}.dose.vcf.gz",
        info = "input/{sample}/chr{chr}.info.gz"
    params:
        pwd = lambda wildcards: SAMPLES.loc[wildcards.sample]['JOB']['pwd'],
        dir = "HRC_impute/{sample}"
    shell: "unzip -P {params.pwd} {input} -d {params.dir}"

rule filter:
    input: rules.unzip.output.vcf
    output: temp("temp/{sample}/chr{chr}.dose.vcf.gz")
    shell:
        """
        ml bcftools
        bcftools view -e "INFO/MAF < 0.01 | INFO/R2 < 0.3" {input} -Oz -o {output}
        bcftools index -t {output}
        """

rule concat:
    input: expand("temp/{{sample}}/chr{chr}.dose.vcf.gz", chr = CHR)
    output: "HRC_impute/final/{sample}_chrall.dose.vcf.gz"
    shell:
        """
        ml bcftools
        bcftools concat {input} -Oz -o {output}
        bcftools index -t {output}
        """

rule snp_info:
    input: rules.concat.output
    output: "HRC_impute/final/{sample}_snp_info.txt"
    shell:
        r"""
        ml bcftools miller
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%AF\t%MAF\t%R2\t%IMPUTED\t%TYPED\t%TYPED_ONLY\n" -o {output} {input}
        """
# mlr --tsv --implicit-csv-header label CHROM,POS,REF,ALT,QUAL,FILTER,AF,MAF,R2,IMPUTED,TYPED,TYPED_ONLY > {output}

rule plink:
    input: rules.concat.output
    output: multiext("HRC_impute/final/{sample}_chrall", '.bed', '.bim', '.fam')
    params: out = "HRC_impute/final/{sample}_chrall"
    shell:
        """
        ml plink
        plink --vcf {input} --double-id --keep-allele-order --make-bed --out {params.out}
        """

rule fix_fam:
    input:
        demo = "../20201015/ADNIMERGE.csv",
        fam = "HRC_impute/final/{sample}_chrall.fam"
    output: out = "HRC_impute/final/{sample}_chrall_fixed.fam"
    script: 'workflow/scripts/fix_fam.R'
