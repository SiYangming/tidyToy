library(readr)
library(dplyr)
riskScore <- read_tsv("data-raw/unIndepInput.txt") %>%
  dplyr::select(id, riskScore)
save(riskScore, file = "data/riskScore.rda")
