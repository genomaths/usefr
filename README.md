<!-- README.md is generated from README.Rmd. Please edit that file -->
USEF-R [<img src="man/figures/logo.png" align="right" />](https://genomaths.github.io/usefr)
==========================================================

<br>

# Welcome to the R package _usefr_

This is an utility package where some useful functions frequently needed to
identify the best fitted probability distribution model for a given data set
are provided.

------------

## Status

This application is under development. Watch this repo or check for updates.

------------

## Dependences

This package depends, so far, from: _BiocParallel_, _minpack.lm_, _numDeriv_,
_copula_. There are also other dependencies which are included in R by default,
e.g., _start_

------------

### Install R dependencies:

```install
    if (!requireNamespace("BiocManager")) install.packages("BiocManager")
    BiocManager::install()
    
    BiocManager::install(c("BiocParallel","minpack.lm", "numDeriv", "copula", 
                            "mclust", "nls2", "cubature", "mixdist"), 
                        dependencies=TRUE)
```

------------

## You can install _*usefr*_ package from GitHub

```install.p
   BiocManager::install("genomaths/usefr")
```


------------

## Tutorials:

<a href="https://genomaths.github.io/usefr/articles/usefr.html/">Get started-with ‘usefr’</a>


<a href="https://genomaths.github.io/usefr/articles/fisher-exact-test-failure-can-lead-to-biassed-results.html/">Fisher’s exact test failure can lead to biased results</a>


------------



