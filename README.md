# R code for the manuscript "Nonparametric quantile regression captures regional variability and scaling deviations in Atlantic surfclam length–weight relationships" 


Manuscript authors: Gorka Bidegain, Marta Sestelo, Patricia L. Luque, Eric N. Powell, Arantza Irirarte, Ibon Uriarte, and Daphne Munroe.

This code was written by Marta Sestelo. In case of questions or comments, please contact sestelo@uvigo.gal. 

This GitHub repository contains the routines required to reproduce the models and testing procedures of the paper entitled “Nonparametric quantile regression captures regional variability and scaling deviations in Atlantic surfclam length–weight relationships”.

To run the cited routines, it is needed to install the npregfast package from [CRAN](https://cran.r-project.org/web/packages/npregfast).

Information about the R Session:
```
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Qtools_1.5.9    np_0.60-17      ggplot2_3.5.2   mgcv_1.9-1     
 [5] nlme_3.1-164    sfsmisc_1.1-17  quantreg_5.97   SparseM_1.81   
 [9] npregfast_1.5.2 patchwork_1.3.1

loaded via a namespace (and not attached):
 [1] sandwich_3.1-0      generics_0.1.4      gtools_3.9.5       
 [4] glmx_0.2-0          lattice_0.22-5      cubature_2.1.0     
 [7] magrittr_2.0.3      grid_4.3.2          RColorBrewer_1.1-3 
[10] iterators_1.0.14    pkgload_1.3.4       foreach_1.5.2      
[13] doParallel_1.0.17   Matrix_1.6-5        conquer_1.3.3      
[16] Formula_1.2-5       survival_3.5-8      scales_1.4.0       
[19] numDeriv_2016.8-1.1 codetools_0.2-19    cli_3.6.5          
[22] rlang_1.1.6         crayon_1.5.3        splines_4.3.2      
[25] withr_3.0.2         tools_4.3.2         parallel_4.3.2     
[28] MatrixModels_0.5-3  dplyr_1.1.4         boot_1.3-30        
[31] quantdr_1.2.2       vctrs_0.6.5         R6_2.6.1           
[34] matrixStats_1.5.0   zoo_1.8-12          lifecycle_1.0.4    
[37] MASS_7.3-60.0.1     shinyjs_2.1.0       pkgconfig_2.0.3    
[40] pillar_1.11.0       gtable_0.3.6        glue_1.8.0         
[43] wesanderson_0.3.7   Rcpp_1.1.0          tibble_3.3.0       
[46] lmtest_0.9-40       tidyselect_1.2.1    rstudioapi_0.17.1
[49] farver_2.1.2        xtable_1.8-4        labeling_0.4.3     
[52] compiler_4.3.2      quadprog_1.5-8
```
