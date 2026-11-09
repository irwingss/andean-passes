# 00_setup.R - install the R packages of the pipeline (study area, sea domains, downloads) and record the session.
# Run from the repository root: Rscript scripts/00_setup.R
need <- c("digest", "geodata", "sf", "terra", "yaml")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
for (p in need) suppressPackageStartupMessages(library(p, character.only = TRUE))
dir.create("outputs", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
