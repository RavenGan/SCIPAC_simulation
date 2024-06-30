# SCIPAC_simulation

## File descriptions
The code to generate the simulation data is in the `./R` folder with explanations as follows. Scheme I, II, and III are the three simulation schemes in the paper 

Gan, D., Zhu, Y., Lu, X., & Li, J. (2024). SCIPAC: quantitative estimation of cell-phenotype associations. Genome Biology, 25(1), 119.

Roughly speaking, Scheme I is the simplest one, containing three cell types with a binary phenotype, Scheme II is more complex with seven cell types from a binary phenotype, and Scheme III is generates single cell RNA sequencing data with ordinal phenotype, mimicking the differentiation process of cells.

### Simulation scheme I
* Scheme_I_scRNA-seq.R: The code to generate the single cell simulation data for scheme I.
* Scheme_I_bulkRNA-seq.R: The code to generate the bulk RNA-seq simulation data for scheme I.

### Simulation scheme II
* Scheme_II_scRNA-seq.R: The code to generate the single cell simulation data for scheme II.
* Scheme_II_bulkRNA-seq.R: The code to generate the bulk RNA-seq simulation data for scheme II.

### Simulation scheme III
* Scheme_III_scRNA-seq.R: The code to generate the single cell simulation data for scheme III.
* Scheme_III_bulkRNA-seq.R: The code to generate the bulk RNA-seq simulation data for scheme III.