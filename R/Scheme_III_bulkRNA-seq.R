rm(list = ls())
set.seed(7)

library(dplyr)
library(splatter)

# This simulation may be time-consuming.
# Consider using a computation cluster

set.seed(1)
sim.scRNA <- readRDS("./Splatter.simulation3/scRNA-seq/sim.raw.scRNA.rds")
sim.scRNA.meta <- readRDS("./Splatter.simulation3/scRNA-seq/sim.raw.scRNA.meta.rds")

all(sim.scRNA.meta$Cell == colnames(sim.scRNA)) # TRUE

# Stage 1: Step >= 0 & Step <= 500
# Stage 2: Step >= 501 & Step <= 1000
# Stage 3: Step >= 1001 & Step <= 1500
# Stage 4: Step >= 1501 & Step <= 2000

# Create four stages
# Stage 1 meta
stage1.idx <- which(sim.scRNA.meta$Step >= 0 & sim.scRNA.meta$Step <= 500)
stage1.meta <- sim.scRNA.meta[stage1.idx, ]
stage1.meta$stages <- "Stage I"

# Stage 2 meta
stage2.idx <- which(sim.scRNA.meta$Step >= 501 & sim.scRNA.meta$Step <= 1000)
stage2.meta <- sim.scRNA.meta[stage2.idx, ]
stage2.meta$stages <- "Stage II"

# Stage 3 meta
stage3.idx <- which(sim.scRNA.meta$Step >= 1001 & sim.scRNA.meta$Step <= 1500)
stage3.meta <- sim.scRNA.meta[stage3.idx, ]
stage3.meta$stages <- "Stage III"

# Stage 4 meta
stage4.idx <- which(sim.scRNA.meta$Step >= 1501 & sim.scRNA.meta$Step <= 2000)
stage4.meta <- sim.scRNA.meta[stage4.idx, ]
stage4.meta$stages <- "Stage IV"


# Sample Stage bulk

sample.bulk <- function(sim.scRNA, stage.meta, n.cases, stage.name){
  save.path <- paste0("./Splatter.simulation3/bulkRNA-seq/", stage.name, ".cases.rds")
  
  stage.scRNA <- sim.scRNA[, stage.meta$Cell]
  stage.col.idx <- c(1:ncol(stage.scRNA))
  stage.cases <- data.frame(matrix(NA, nrow = nrow(stage.scRNA), ncol = n.cases))
  
  for (i in 1:n.cases) {
    set.seed(i)
    cell.idx <- sample(stage.col.idx, 10000, replace = TRUE)
    sample.i <- stage.scRNA[, cell.idx]
    stage.cases[, i] <- rowSums(sample.i)/ncol(sample.i)
  }
  rownames(stage.cases) <- rownames(stage.scRNA)
  saveRDS(stage.cases, save.path)
}

sample.stageI <- sample.bulk(sim.scRNA, stage1.meta, 50, "StageI")
sample.stageII <- sample.bulk(sim.scRNA, stage2.meta, 50, "StageII")
sample.stageIII <- sample.bulk(sim.scRNA, stage3.meta, 50, "StageIII")
sample.stageIV <- sample.bulk(sim.scRNA, stage4.meta, 50, "StageIV")







