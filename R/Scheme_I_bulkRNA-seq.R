rm(list = ls())
set.seed(7)

library(dplyr)
library(splatter)

set.seed(1)
sim.scRNA <- readRDS("./Splatter.simulation/scRNA-seq/sim.raw.scRNA.rds")
sim.scRNA.meta <- readRDS("./Splatter.simulation/scRNA-seq/sim.raw.scRNA.meta.rds")

# Create two phenotypes.
group1.idx <- which(sim.scRNA.meta$Group == "Group1")
group1.pheno1.idx <- sample(group1.idx, size = floor(length(group1.idx)/2))
group1.pheno2.idx <- group1.idx[group1.idx %in% group1.pheno1.idx == FALSE]

# Phenotype 1
group1.pheno1 <- sim.scRNA.meta[group1.pheno1.idx, ]
group2.pheno1.idx <- which(sim.scRNA.meta$Group == "Group2")
group2.pheno1 <- sim.scRNA.meta[group2.pheno1.idx, ]

pheno1.meta <- rbind(group1.pheno1, group2.pheno1)
pheno1.meta$Phenotype <- "Phenotype.1"

# Phenotype2
group1.pheno2 <- sim.scRNA.meta[group1.pheno2.idx, ]
group3.pheno2.idx <- which(sim.scRNA.meta$Group == "Group3")
group3.pheno2 <- sim.scRNA.meta[group3.pheno2.idx, ]

pheno2.meta <- rbind(group1.pheno2, group3.pheno2)
pheno2.meta$Phenotype <- "Phenotype.2"

# Sample phenotype 1 bulk
pheno1.scRNA <- sim.scRNA[, pheno1.meta$Cell]
pheno1.col.idx <- c(1:ncol(pheno1.scRNA))
n.cases <- 50
pheno1.cases <- data.frame(matrix(NA, nrow = nrow(pheno1.scRNA), ncol = n.cases))

for (i in 1:n.cases) {
  set.seed(i)
  cell.idx <- sample(pheno1.col.idx, 10000, replace = TRUE)
  sample.i <- pheno1.scRNA[, cell.idx]
  pheno1.cases[, i] <- rowSums(sample.i)/ncol(sample.i)
}
rownames(pheno1.cases) <- rownames(pheno1.scRNA)
saveRDS(pheno1.cases, "./Splatter.simulation/bulkRNA-seq/pheno1.cases.rds")

# Sample phenotype 2 bulk
pheno2.scRNA <- sim.scRNA[, pheno2.meta$Cell]
pheno2.col.idx <- c(1:ncol(pheno2.scRNA))
n.cases <- 50

pheno2.cases <- data.frame(matrix(NA, nrow = nrow(pheno2.scRNA), ncol = n.cases))

for (i in 1:n.cases) {
  set.seed(i)
  cell.idx <- sample(pheno2.col.idx, 10000, replace = TRUE)
  sample.i <- pheno2.scRNA[, cell.idx]
  pheno2.cases[, i] <- rowSums(sample.i)/ncol(sample.i)
}
rownames(pheno2.cases) <- rownames(pheno2.scRNA)
saveRDS(pheno2.cases, "./Splatter.simulation/bulkRNA-seq/pheno2.cases.rds")
