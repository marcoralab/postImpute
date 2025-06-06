#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
library(vroom)
library(tidyr)
library(purrr)

if (!exists("snakemake")) {
  setClass("snakemake_fake", representation(
    params = "list", input = "list", output = "list",
    log = "list", wildcards = "list"))
  snakemake <- new("snakemake_fake",
    input = list(
      oldfam = c("temp/raw/GSA1-AD1_single-probe.fam",
                 "temp/raw/GSA1-AD2_single-probe.fam",
                 "temp/raw/GSA1-PD1_single-probe.fam",
                 "temp/raw/GSA1-PD2_single-probe.fam"),
      newfam = "temp/raw/raw.fam"
    ),
    params = list(),
    log = list(),
    output = list("results/raw/raw.fam"),
    wildcards = list()
  )
}

col.n <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.nn <- c("none", "FIDIID", "PID_new", "MID_new", "Sex_new", "Phe_new")
col.t <- "ccccii"

new_fam <- snakemake@input[["newfam"]] %>%
  vroom::vroom(col_names = col.nn, col_types = col.t)

stopifnot(nrow(distinct(new_fam, FIDIID)) == nrow(new_fam))

old_fam <- snakemake@input[["oldfam"]] %>%
  map_dfr(vroom::vroom, col_names = col.n, col_types = col.t) %>%
  unite("FIDIID", FID, IID, sep = "_", remove = F) %>%
  distinct(FIDIID, .keep_all = TRUE)

stopifnot(nrow(distinct(old_fam, FIDIID)) == nrow(old_fam))

left_join(new_fam, old_fam, by = "FIDIID") %>%
  mutate(FID = ifelse(is.na(FID), FIDIID, FID)) %>%
  mutate(IID = ifelse(is.na(IID), FIDIID, IID)) %>%
  mutate(PID = ifelse(is.na(PID), PID_new, PID)) %>%
  mutate(MID = ifelse(is.na(MID), MID_new, MID)) %>%
  mutate(Sex = ifelse(is.na(Sex), Sex_new, Sex)) %>%
  mutate(Phe = ifelse(is.na(Phe), Phe_new, Phe)) %>%
  select(!!col.n) %>%
  mutate(IID = gsub("�", "", IID)) %>%
  vroom::vroom_write(snakemake@output[[1]], col_names = F, delim = "\t",
                     quote = "none")
