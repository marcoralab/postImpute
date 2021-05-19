#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(magrittr)

fname <- snakemake@input[['vcf']]
mapfile <- snakemake@input[['mapping']]
foutmap <- snakemake@output[['mapping']]
foutrehead <- snakemake@output[['rehader']]

conv_sampline <- function (sampline) {
  sampline %>%
    strsplit("\t") %>%
    unlist %>%
    tibble(ID = .) %>%
    slice(10:n())
}

read_samples_raw <- function (vcf) {
  while (!exists("samples")) {
    line <- readLines(vcf, 1)
    if (grepl("^#CHROM", line)) {
      samples <- conv_sampline(line)
    }
  }
  return(samples)
}

read_samples_preproc <- function (sl) {
  sampline <- readLines(sl, 1)
  if (grepl("^#CHROM", sampline)) {
    samples <- conv_sampline(sampline)
  } else {
    stop("Sample line not found.")
  }
  return(samples)
}

read_samples <- function (fname) {
  sampline <- "zcat %s | head -n200 | grep \"^#CHROM\"" %>%
    sprintf(fname) %>%
    system(intern = T)
  if (grepl("^#CHROM", sampline)) {
    samples <- conv_sampline(sampline)
  } else {
    stop("Sample line not found.")
  }
  return(samples)
}

samples <- read_samples(fname)

mapping <- mapfile %>%
  read_tsv(col_types = cols(.default = "c"))

find.col.same <- function(famcol) {
  samplist <- pull(samples, famcol)
  countsame <- function (mapcol) {
    sum(samplist %in% pull(mapping, mapcol))
  }
  cts <- sapply(names(mapping), countsame)
  max_cts <- names(cts[cts == max(cts)])
  if ( max(cts) < dim(samples)[1]) {
    print(cts)
    best <- max_cts[1]
    mismaps <- samplist[!(samplist %in% pull(mapping, best))]
    print(mismaps)
    stop("Cannot map all IDs")
  }
  max_cts[1]
}

mapping %<>%
  unite(raw, raw_FID, raw_IID, sep = "_") %>%
  unite(old, old_FID, old_IID, sep = "_") %>%
  unite(oldest, oldest_FID, oldest_IID, sep = "_")

bycol <- find.col.same("ID")
mapped <- samples %>%
  rename(!!bycol := ID) %>%
  left_join(mapping, by = bycol) %>%
  mutate(FID_IID = paste0(FID, "_", IID)) %>%
  rename(ID = !!bycol) %>%
  select(FID, IID, FID_IID, ID)

stopifnot(!any(is.na(mapped$FID) | is.na(mapped$IID)))
stopifnot(!any(grepl("_", mapped$FID) | grepl("_", mapped$IID)))

write_tsv(mapped, path = foutmap)
mapped %>%
  select(FID_IID) %>%
  write_tsv(path = foutrehead, col_names = F)

message("Finished mapping sample names")
