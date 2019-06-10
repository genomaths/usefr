## Welcome to the R package _usefr_

This is an utility package where some useful functions frequently needed in the dowstream statistical analyses will be available.

------------
### Status

   This application is under development. Watch this repo or check for updates.

------------
### Dependences

This package depends, so far, from: _BiocParallel_, _minpack.lm_, _numDeriv_, _copula_. There are also other dependencies which are included in R by default, e.g., _start_


------------
### Install R dependencies:

```install
    if (!requireNamespace("BiocManager")) install.packages("BiocManager")
    BiocManager::install()
    
    BiocManager::install('BiocParallel')
    install.packages(c("minpack.lm", "numDeriv", "copula"), dependencies=TRUE)
```

------------
### You can install _*usefr*_ package from GitHub

```install.p
   devtools::install_git("https://github.com/genomaths/usefr.git")

```

------------
### _usefr_ R Package Manual:

<a href="https://github.com/genomaths/usefr/blob/master/usefr.pdf" target="_blank">_usefr_ PDF manual</a>


<a href="https://genomaths.github.io/usefr_manual/usefr_manual.html" target="_blank">_usefr_ browser manual</a>

