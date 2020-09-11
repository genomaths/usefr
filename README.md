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
the package source. 

"/lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by ...
cubature.so"

Or

"./src/divonne/Split.c:119:3: note: use option -std=c99 or -std=gnu99 to compile your code"

In the above situations, just proceed with the following steps (<https://github.com/bnaras/cubature/issues/29>):
  
  1. In your home directory search for the folder named ".R" (a hidden folder) and a file named "Makevars", i.e.: "~/.R/Makevars".
     If you do not have the folder ".R", then creates it. If you do not have the file "Makevars", 
     then you can create an empty text file 
     
  2. Adds the following line to the the file "Makevars":
  
     ```
          CFLAGS=-std=gnu99
     ```
  
     if your file is empty, the you can add, e.g., something like this:
  
     ```
     SOLARIS=$(shell $(R_HOME)/bin/Rscript -e 'cat(grepl("SunOS", Sys.info()["sysname"]))')
     ifeq ($(SOLARIS),TRUE)
     SOLARIS_FLAG=-DSOLARIS
     else
     SOLARIS_FLAG=-USOLARIS
     endif
     
     CFLAGS=-std=gnu99
     
     ```

If the above step does not works then,  follow the link: <https://github.com/cdr/code-server/issues/347>

------------

### You can install _*usefr*_ package from GitHub

```install.p
   devtools::install_git("https://github.com/genomaths/usefr.git")

```

Or download the binary of 'usefr' R package from here: <https://is.gd/4BKQQN>
and in the R console type:

```install.usefr

install.packages("usefr_0.1.0_R_x86_64-pc-linux-gnu.tar.gz", 
                 repos = NULL, type = "source")

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



