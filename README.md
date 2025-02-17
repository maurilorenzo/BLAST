# BLAST

This repository contains the Supporting R/Cpp code for the article "Spectral decomposition-assisted multi-study factor analysis".

The repo is structured as follows:

```
├── application
    ├── data
        ├── genedata.rds
    ├── gene_preprocessing.R
    ├── genes_model_comparison.Rmd
    ├── genes_plot.R
    └── helpers_application.R
├── simulations
    ├── helpers_numerical_experiments.R
    ├── numerical_experiments.R
    └── numerical_experiments_supplemental.R
├── blast_wrapper.R
├── example_vignette.Rmd
└── helper_functions.cpp
 ```  

The `blast_wrapper.R` script contain an R wrapper for implementing the methodology. `helper_functions.cpp` contains the Cpp subroutines for the method.  
The `simulations` folder contains the code to replicate the numerical experiments. In particular, `numerical_experiments.R` and `numerical_experiments_supplemental.R` contain code to replicate the numerical experiments presented in Section 4 and Supplementary Material respectively.  
The `application` folder contains the code to replicate the application on gene expression data. The notebook `genes_model_comparison.Rmd` replicates the results presented in Section 5 and the script `genes_plot.R` produces the plots presented in the Supplementary Material.  
`example_vignette.Rmd` presents a simple example of our methodology



