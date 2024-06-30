rm(list = ls())
set.seed(7)

library(dplyr)
library(splatter)
library(scater)
set.seed(1)
sim1 <- splatSimulate(group.prob = c(0.25, 0.25, 0.25, 0.25), method = "paths",
                      de.prob = 0.5,
                      path.from = c(0, 0, 0, 0),
                      verbose = FALSE,
                      batchCells = 10000,
                      path.nSteps = 2000)

sim2 <- logNormCounts(sim1)
sim2 <- runPCA(sim2)
pdf("./Splatter.simulation3/res.Oct.04.2022/sim.plot.group.pdf", height = 7, width = 7)
plotPCA(sim2, colour_by = "Group") + ggtitle("Four groups")
dev.off()

pdf("./Splatter.simulation3/res.Oct.04.2022/sim.plot.step.pdf", height = 7, width = 7)
plotPCA(sim2, colour_by = "Step") + ggtitle("Four groups")
dev.off()

exprs.dat <- sim1@assays@data@listData[["counts"]]
meta.data <- sim1@colData@listData %>% as.data.frame()

saveRDS(exprs.dat, "./Splatter.simulation3/scRNA-seq/sim.raw.scRNA.rds")
saveRDS(meta.data, "./Splatter.simulation3/scRNA-seq/sim.raw.scRNA.meta.rds")
