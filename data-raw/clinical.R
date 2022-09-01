library(readr)
library(dplyr)
clinical <- read_tsv("data-raw/unIndepInput.txt")
save(clinical, file = "data/clinical.rda")

