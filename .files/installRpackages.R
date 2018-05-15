packages <- c("tidyverse", "magrittr", "stringr",
              "ggplot2", "scales", "feather",
              "kableExtra", "doParallel")
for (i in packages) {
  missing <- suppressWarnings(!require(i,
    character.only = TRUE, quietly = TRUE))
  if (missing) {
    install.packages(i, repos = "https://mirrors.sorengard.com/cran/")
  }
}
