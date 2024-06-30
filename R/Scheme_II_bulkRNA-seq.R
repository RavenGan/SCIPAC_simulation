rm(list = ls())
set.seed(7)

library(dplyr)
library(splatter)

set.seed(1)
sim.scRNA <- readRDS("./Splatter.simulation2/scRNA-seq/sim.raw.scRNA.rds")
sim.scRNA.meta <- readRDS("./Splatter.simulation2/scRNA-seq/sim.raw.scRNA.meta.rds")

# Phenotype 1: Group1, Group3, half of Group5, 6, 7
# Phenotype 2: Group2, Group4, half of Group5, 6, 7


# Create two phenotypes.
group5.idx <- which(sim.scRNA.meta$Group == "Group5")
group5.pheno1.idx <- sample(group5.idx, size = floor(length(group5.idx)/2))
group5.pheno2.idx <- group5.idx[group5.idx %in% group5.pheno1.idx == FALSE]

group6.idx <- which(sim.scRNA.meta$Group == "Group6")
group6.pheno1.idx <- sample(group6.idx, size = floor(length(group6.idx)/2))
group6.pheno2.idx <- group6.idx[group6.idx %in% group6.pheno1.idx == FALSE]

group7.idx <- which(sim.scRNA.meta$Group == "Group7")
group7.pheno1.idx <- sample(group7.idx, size = floor(length(group7.idx)/2))
group7.pheno2.idx <- group7.idx[group7.idx %in% group7.pheno1.idx == FALSE]

# Phenotype 1
group567.pheno1 <- sim.scRNA.meta[c(group5.pheno1.idx, 
                                    group6.pheno1.idx,
                                    group7.pheno1.idx), ]
group1.pheno1.idx <- which(sim.scRNA.meta$Group == "Group1")
group1.pheno1 <- sim.scRNA.meta[group1.pheno1.idx, ]

group3.pheno1.idx <- which(sim.scRNA.meta$Group == "Group3")
group3.pheno1 <- sim.scRNA.meta[group3.pheno1.idx, ]

pheno1.meta <- rbind(group1.pheno1, group3.pheno1) %>%
  rbind(group567.pheno1)
pheno1.meta$Phenotype <- "Phenotype.1"

# Phenotype 2
group567.pheno2 <- sim.scRNA.meta[c(group5.pheno2.idx, 
                                    group6.pheno2.idx,
                                    group7.pheno2.idx), ]
group2.pheno2.idx <- which(sim.scRNA.meta$Group == "Group2")
group2.pheno2 <- sim.scRNA.meta[group2.pheno2.idx, ]

group4.pheno2.idx <- which(sim.scRNA.meta$Group == "Group4")
group4.pheno2 <- sim.scRNA.meta[group4.pheno2.idx, ]

pheno2.meta <- rbind(group2.pheno2, group4.pheno2) %>%
  rbind(group567.pheno2)
pheno2.meta$Phenotype <- "Phenotype.2"


# Sample phenotype 1 bulk
pheno1.scRNA <- sim.scRNA[, pheno1.meta$Cell]
pheno1.col.idx <- c(1:ncol(pheno1.scRNA))
n.cases <- 100
pheno1.cases <- data.frame(matrix(NA, nrow = nrow(pheno1.scRNA), ncol = n.cases))

for (i in 1:n.cases) {
  set.seed(i)
  cell.idx <- sample(pheno1.col.idx, 10000, replace = TRUE)
  sample.i <- pheno1.scRNA[, cell.idx]
  pheno1.cases[, i] <- rowSums(sample.i)/ncol(sample.i)
}
rownames(pheno1.cases) <- rownames(pheno1.scRNA)
saveRDS(pheno1.cases, "./Splatter.simulation2/bulkRNA-seq/pheno1.cases.rds")

# Sample phenotype 2 bulk
pheno2.scRNA <- sim.scRNA[, pheno2.meta$Cell]
pheno2.col.idx <- c(1:ncol(pheno2.scRNA))
n.cases <- 100

pheno2.cases <- data.frame(matrix(NA, nrow = nrow(pheno2.scRNA), ncol = n.cases))

for (i in 1:n.cases) {
  set.seed(i)
  cell.idx <- sample(pheno2.col.idx, 10000, replace = TRUE)
  sample.i <- pheno2.scRNA[, cell.idx]
  pheno2.cases[, i] <- rowSums(sample.i)/ncol(sample.i)
}
rownames(pheno2.cases) <- rownames(pheno2.scRNA)
saveRDS(pheno2.cases, "./Splatter.simulation2/bulkRNA-seq/pheno2.cases.rds")

