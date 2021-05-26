
rule make_bgen:
    input:
        gen = renamed
    output:
        bgen = temp("temp/{sample}_chr{chrom}_filtered.bgen"),
        samp = temp("temp/{sample}_chr{chrom}.sample")
    shell:
        """
module load qctool/v2
qctool -g {input.gen} -vcf-genotype-field GP \
-os {output.samp} -og {output.bgen}
"""

rule cat_bgen_samp:
    input:
        gen = expand("temp/{{sample}}_chr{chrom}_filtered.bgen", chrom=CHROM),
        samp = expand("temp/{{sample}}_chr{chrom}.sample", chrom=CHROM)[0]
    output:
        gen = "data/{sample}_chrall_filtered.bgen",
        samp = "data/{sample}_chrall.sample"
    shell:
        """
cat-bgen -g {input.gen} -og {output.gen}
cp {input.samp} {output.samp}"""

bga_gen = expand("temp/{sample}_chr{{chrom}}_filtered.bgen", sample=SAMPLE)
bga_samp = expand("temp/{sample}_chr{{chrom}}.sample", sample=SAMPLE)

rule make_bgen_allsamp:
    input:
        gen = bga_gen,
        samp = bga_samp
    output:
        bgen = "data/merged/merged_chr{chrom}_filtered.bgen",
        samp = "data/merged/merged_chr{chrom}.sample"
    params:
        args = " ".join(["-g {} -s {}".format(gen, samp) for gen, samp in zip(bga_gen, bga_samp)])
    threads: 10
    shell:
        """
module load qctool/v2
qctool {params.args} -og {output.bgen} -os {output.samp} -threads 10
"""

rule cat_bgen_allsamp:
    input:
        gen = expand("data/merged/merged_chr{chrom}_filtered.bgen", chrom=CHROM),
        samp = expand("data/merged/merged_chr{chrom}.sample", chrom=CHROM)[0]
    output:
        gen = "data/merged/merged_chrall_filtered.bgen",
        samp = "data/merged/merged_chrall.sample"
    shell: "cat-bgen -g {input.gen} -og {output.gen}; cp {input.samp} {output.samp}"