import pysam

with pysam.VariantFile(snakemake.input[0], 'r') as vcf:
    info = vcf.header.info
    filt = vcf.header.filters

if 'IMPUTED' in info:
    ver = '4'
elif 'GENOTYPED' in filt:
    ver = '3'
else:
    raise  Exception("Could not detect minimac version")

with open(snakemake.output[0], 'w') as f:
    print(ver, file=f)