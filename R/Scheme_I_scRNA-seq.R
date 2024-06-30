rm(list = ls())
set.seed(7)

library(dplyr)
library(splatter)

set.seed(1)
sim1 <- splatSimulate(group.prob = c(0.8, 0.1, 0.1), method = "groups",
                      de.prob = c(0.1, 0.01, 0.01),
                      verbose = FALSE,
                      batchCells = 10000)
exprs.dat <- sim1@assays@data@listData[["counts"]]
meta.data <- sim1@colData@listData %>% as.data.frame()

saveRDS(exprs.dat, "./Splatter.simulation/scRNA-seq/sim.raw.scRNA.rds")
saveRDS(meta.data, "./Splatter.simulation/scRNA-seq/sim.raw.scRNA.meta.rds")
