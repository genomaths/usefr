## Welcome to the R package _usefr_

This is an utility package where some useful functions frequently needed in the dowstream statistical analyses will be available.
The package started with the first two functions and other functions will be added "little by little".

### Status

   This application is under development. Watch this repo or check for updates.

### Dependences

This package depends, so far, from: _BiocParallel_, _minpack.lm_, _numDeriv_. There are also other dependencies which are included in R by default, e.g., _start_


### Install R dependencies:

```install
    source("https://bioconductor.org/biocLite.R")
    biocLite('BiocParallel')
    install.packages(c("minpack.lm", "numDeriv"),dependencies=TRUE)

```

### You can install _*usefr*_ package from GitHub

```install.p
   devtools::install_git("https://github.com/genomaths/usefr.git")

```
