library(tidyverse)
# https://www.cell.com/cms/10.1016/j.immuni.2013.10.003/attachment/8dc04d32-6aff-4eda-99f5-6401bcae9751/mmc1.pdf
immunity <- read_csv("data-raw/immunity/immune_cell_gene.csv")
idx <- !immunity$CellType %in% c("Blood vessels", "Normal mucosa", "SW480 cancer cells", "Lymph vessels")
immunity <- immunity[idx, ]
immunity <- immunity %>%
  split(., .$CellType) %>%
  lapply(., function(x) (x$Symbol))
immunity24 <- lapply(immunity, unique)
save(immunity24, file = "data/immunity24.rda")
write_csv(as.data.frame(names(immunity24)),
          file = "data-raw/immunity/immunity24cells.csv",
          col_names = FALSE)
