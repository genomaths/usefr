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
    install.packages(c("minpack.lm", "numDeriv", "copula", "mclust",
                        "nls2", "cubature"), dependencies=TRUE)
```

### Installing issues with the R package "cubature" on CentOS

The installation of "cubature" R package on CenOS (required by 'usefr') from 
the source, would produce some error when trying to compile the C++ code from
the package source. Use instead the binary package, which is already compiled.
After that, still some error message can be returned, like this: 

"/lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by ...
cubature.so"

In the above situation, just check one of the following alternatives:
  
  1. Download the binary of 'cubature' R package from here: <https://is.gd/WZI32A> 
     and then the binary of 'usefr' R package from here: <https://is.gd/4BKQQN>.
     
  2. Or follow the link: <https://github.com/cdr/code-server/issues/347>  
  
  3. Or just install [Anaconda](https://docs.anaconda.com/anaconda/install/linux/) in your CentOS and then [install R on Annaconda](https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/) and all the nightmares with R on CentOS will pass away!

------------

### You can install _*usefr*_ package from GitHub

```install.p
   devtools::install_git("https://github.com/genomaths/usefr.git")

```

------------
### _usefr_ R Package Manual:

<a href="https://github.com/genomaths/usefr/blob/master/usefr.pdf" target="_blank">_usefr_ PDF manual</a>


<a href="https://genomaths.github.io/usefr_manual/usefr_manual.html" target="_blank">_usefr_ browser manual</a>

------------

### Some examples of applications:
<a href="https://genomaths.com/stats/sampling-from-a-mixture-of-distributions/">Sampling from a Mixture of Distributions</a>


<a href="https://genomaths.com/stats/non-linear-fit-of-mixture-distributions/">Fitting Mixture Distributions</a>


------------



