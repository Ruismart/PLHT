# PLHT

## Description     

A ***P***CA-based ***L***inear ***H***olistic ***T***rajectory method to measure the postnatal development level for mouse intestinal Epithelium/Fibroblast: 

Using mouse intestinal scRNAseq WT data (Epcam+ Epithelial & Pdgfra+ Mesenchymal) and trajectory features merged/filtered from stage-DEGs (P00-P07-P14-P21-P28), we have built a PCA-based linear model to regress A holistical time coordinate along postnatal development, which could be applied to the comparison of maturity level between CTL-vs-CKO, as well within any other normalized singlecell or bulk transcriptome datasets.


## Installation       

You can install the initial version of PLHT like so:

``` r
devtools::install_github("Ruismart/PLHT", subdir="v1.0.0")
```

## Tutorial            

### This is a basic example:

``` r
# load
library(PLHT)

## basic example code
# for bulk RNAseq, input matrix
traj.test  <- plht(mat = log2(mat_pc+1),
                   design = design,
                   stsp = "KO",       # set standard sample condition
                   stsp.type = "median",     # would normalize median of standard to base level
                   scale = TRUE,
                   mat.coef.EpC = mat.coef.EpC,
                   mat.coef.Fib = mat.coef.Fib,
                   features.overlap = T)


# for scRNAseq, input seurat obj
test.seur <- plht.seur(seur.obj = test.seur,
                       assay = "RNA",
                       condition = "condition",
                       stsp = "P0",
                       stsp.type = "median",
                       mat.coef.EpC = mat.coef.EpC,
                       mat.coef.Fib = mat.coef.Fib
                       )
```

### More details see links below: 


#### Report of modeling
  &emsp;&emsp;[**PLHT.modeling**](https://github.com/Ruismart/PLHT/tree/main/tutorial/)       


#### Tutorial using several public datasets           

  &emsp;&emsp;[**local bulk**](https://github.com/Ruismart/PLHT/tree/main/tutorial/) 

  &emsp;&emsp;[**public single**](https://github.com/Ruismart/PLHT/tree/main/tutorial/) 



## 

--

## How to cite       

...wait...











