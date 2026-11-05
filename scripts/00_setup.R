# 00_setup.R - install the R packages of the pipeline and record the session.
# Run from the repository root: Rscript scripts/00_setup.R
need <- c("ape", "data.table", "digest", "geodata", "geosphere", "ggplot2", "ggrepel", "jsonlite", "maps", "mgcv",
          "openxlsx", "patchwork", "png", "ragg", "rayshader", "rgl", "scales", "sf", "terra", "yaml")
miss <- setdiff(need, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
for (p in need) suppressPackageStartupMessages(library(p, character.only = TRUE))
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
