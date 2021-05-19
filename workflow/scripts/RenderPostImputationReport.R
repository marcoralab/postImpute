
message("Render Final Report",
        "\n markdown: ", snakemake@input[["markdown"]],
        "\n intermediates_dir: ", snakemake@params[["output_dir"]],
        "\n output_file: ", snakemake@output[["outfile"]],
        "\n output_dir: ", snakemake@params[["output_dir"]],
        "\n path: ", snakemake@params[["path"]],
        "\n chrom: ", snakemake@params[["chrom"]],
        "\n cohort: ", snakemake@params[["cohort"]],
        "\n maf: ", snakemake@params[["maf"]],
        "\n rsq: ", snakemake@params[["rsq"]],
        "\n rsq2: ", snakemake@params[["rsq2"]],
        "\n sampsize: ", snakemake@params[["sampsize"]],
        "\n out: ", snakemake@params[["out"]],
        "\n\n"
      )

rmarkdown::render(
  input = snakemake@input[["markdown"]],
  clean = TRUE,
  intermediates_dir = snakemake@params[["output_dir"]],
  output_file = snakemake@output[["outfile"]],
  output_dir = snakemake@params[["output_dir"]],
  output_format = "all",
  params = list(
    rwd = snakemake@params[['rwd']],
    path = snakemake@params[["path"]],
    chrom = snakemake@params[["chrom"]],
    cohort = snakemake@params[["cohort"]],
    maf = snakemake@params[["maf"]],
    rsq = snakemake@params[["rsq"]],
    rsq2 = snakemake@params[["rsq2"]],
    sampsize = snakemake@params[["sampsize"]],
    # out = snakemake@params[["out"]],
    outpath = snakemake@params[["output_dir"]]
  )
)
