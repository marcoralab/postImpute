# imputePrep

A Snakemake pipeline to process the output from the Michigan Imputation Server.

## Requirements:

 * Anaconda or Miniconda python
   * (Alternatively make your own Python2 environment and copy or symlink to the `scripts/checkEnv` folder of your install location. You will also have to manually copy the other files.)
 * Python3 (The future is now!)
 * Snakemake (Please install first.)

## Installation:

You will need the `install` script and the hidden `.files` folder.

 * Run `install` to install to your present folder.
 * Run `install -p [install prefix]` to install in a different folder.
 * Use the `-d` flag to also delete the install files.

## Usage:

### Config files:

There are two config files with settings for imputation prep:

 * *config.yaml* contains settings for the QC and processing.
 * *cluster.yaml* contains settings for cluster execution.

Please review those files before starting.

### Data-file location:

By default, the pipeline looks for PLINK filesets in the base directory of your installation. You can choose a different destination by editing `config.yaml`.

### Running:

You can excecute the pipeline on a LSF cluster by navigating to the install folder and running `./run -j [number of desired threads] [additional options]`. You can alternatively run the pipeline with `snakemake` and whatever options you wish to use.
