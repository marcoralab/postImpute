# postImpute

A Snakemake pipeline to process the output from the Michigan/TOPMed Imputation Servers.

## Requirements:

* Anaconda or Miniconda python
* Singularity
* Python3 (The future is now!)
* Snakemake

## Usage:

### Configuration:

`config/config.yaml` contains settings for post-imputation processing.

Please review and edit the configuration file before starting.

#### Data-file location:

By default, the pipeline looks for directories containing either unzipped imputation server outputs or zip files in the base directory of your installation. You can choose a different destination by editing `directory_in`. Cohort names are infered from the stem of the plink filesets.

If you do not want all filesets in the `directory_in` to be processed and imputed and the downloads have already been unzipped, then you can include `COHORT:` or `SAMPLES:` and a list of cohorts in the config file like this:

```yaml
COHORT:
  - group1
  - study2
```

If you have not unzipped the downloads, provide a hash, with cohort names as the keys, and job id (`id`) and password (`pwd`) as nested hashes. An example is below:

```yaml
COHORT:
  COHORT_NAME:
    id: "job-20210518-083345-127"
    pwd: "5vZMb8lF\\|TdqSH"
```

You MUST provide this if you have not unzipped your files.

#### Output directory:

The output files for the pipeline will be stored in the directory specified by `out_dir`. This directory *must* exist.

#### Output file choices

Output file choices can be specified under `outputs`, and you can delete or comment out any list item you do not want. The options follow:

  - `stat_report` A report of imputation quality
  - `vcf_bycohort` Bgzipped VCF files for each cohort
  - `vcf_merged` A bgzipped VCF with all cohorts
  - `bgen_bycohort` Binary oxford filesets for each cohort
  - `bgen_merged` A binary oxford fileset with all cohorts
  - `plink_bycohort` Binary plink filesets for each cohort
  - `plink_merged` A binary plink fileset with all cohorts

#### Chromosome selection

Select the chromosomes you want to upload by editing `chroms`. Separate ranges with ":" and listed chromosomes with ",". You can put both in the same string. Use M for mitochondra.

Options are 1:22,X,Y,M

#### Reference and genome build

The pipeline will reference-align the input files before imputation using the fasta file specified under `ref:`. This fasta MUST be in the same genome build as the input. Input genome builds must match the builds supported by the chosen server, and you must specify the build under `imputation` if it does not match the default for the server.

This pipeline supports all builds for the downloaded imputed files. Be aware that TOPMed imputation returns GRCh38 genotypes.

#### QC

The pipeline can filter the files provided by the servers on a per-cohort basis. The options are under `qc:` and are documented in the configuration file.

#### Sample inclusion or excusion by name

You can also include OR exclude samples by a named list file using one of the following options:

- `include_samp:` filename of samples to include
- `exclude_samp:` filename of samples to exclude

`include_samp` accepts a space/tab-delimited text file with vcf IDs in the first column and optional sex in the second column, and removes all unlisted samples from the current analysis. `exclude_samp` does the same for all *listed* samples.

### Running

You must run the pipeline with Anaconda environments and Singularity enabled. You can either use a Snakemake profile to do so or run the pipeline with the following command:

```bash
snakemake --use-conda --use-singularity
```

Make sure your cluster supports the amount of memory required for report generation (528 GiB) and that the memory is being properly requested. If that is not the case, you can edit the resource requirements on the rule `stats`. We have included a `lsf.yaml` file that ensures those resources are available within this pipeline.
