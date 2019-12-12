# Using packrat

[Packrat](https://cran.r-project.org/web/packages/packrat/index.html) allows us to build isolated, reproducible and portable libraries of packages used per project.
It employs a repository-like paradigm where, any instance of R created within the project
has access to only the packages installed in the project's packrat repository.  The packrat repository
can be moved to new locations/platforms and restored.

## git

packrat will create a subdirectory of the project, called `packrat`.  Eventually, the git repository
(different than the packrat repository) will need to include the packrat subdirectory.  But for 
now (during development) we can exclude it.


## packrat initial setup

We found success by initializing an empty packrat, and then manually adding packages.

Load the module with the version of R you want

```
$ module use /mod/bigelow
$ module load R/dada2      #loads R 3.5.1     or is is dada2/1.12 which loads R 3.6 ??
```

Set your working directory to the packrat project and invoke R (without ``--vanilla` I have learned)

```
$ cd /path/to/my/project
$ R
```

Once R is fired up, set package repositories initialize a packrat repos but prevent auto installation of packages

```
options(repos=c(getOption("repos"), "https://cran.mtu.edu/", BiocManager::repositories()))
# prevent packrat::init form scanning any source code
packrat::init(infer.dependencies = FALSE)
# it will only have packrat, but take a snapshot anyway
packrat::snapshot(infer.dependencies = FALSE)
```

Now continue with other packages - note the snapshot partway through, it's not required
but what the  heck.

```
install.packages("dplyr")
install.packages("readr")
install.packages("configr")
install.packages("futile.logger")
install.packages("BiocManager")
# update the snapshot
packrat::snapshot(infer.dependencies = FALSE)

install.packages("Biostrings")
install.packages("ShortRead")  
install.packages("dada2")
# update the snapshot
packrat::snapshot(infer.dependencies = FALSE)
```


## packrat - a peek inside

Let's see what we have in our project library.  Note they are within the project packrat library.


```
x <- dplyr::as_tibble(installed.packages())
print(x, %>% dplyr::select(Package, LibPath, Version), n = nrow(x))
# # A tibble: 102 x 3
#     Package          LibPath                                            Version 
#     <chr>            <chr>                                              <chr>   
#   1 assertthat       /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.2.1   
#   2 backports        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.1.5   
#   3 BH               /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.69.0-1
#   4 Biobase          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.42.0  
#   5 BiocGenerics     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.28.0  
#   6 BiocManager      /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.30.10 
#   7 BiocParallel     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.16.6  
#   8 Biostrings       /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.50.2  
#   9 bitops           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.0-6   
#  10 cli              /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.0.0   
#  11 clipr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.7.0   
#  12 colorspace       /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4-1   
#  13 configr          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.4   
#  14 crayon           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.3.4   
#  15 dada2            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.10.1  
#  16 data.table       /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.12.8  
#  17 DelayedArray     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.8.0   
#  18 digest           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.6.23  
#  19 dplyr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.8.3   
#  20 ellipsis         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.0   
#  21 fansi            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.4.0   
#  22 farver           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.0.1   
#  23 formatR          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.7     
#  24 futile.logger    /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4.3   
#  25 futile.options   /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.0.1   
#  26 GenomeInfoDb     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.18.2  
#  27 GenomeInfoDbData /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.2.0   
#  28 GenomicAlignmen… /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.18.1  
#  29 GenomicRanges    /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.34.0  
#  30 ggplot2          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 3.2.1   
#  31 glue             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.3.1   
#  32 gtable           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.0   
#  33 hms              /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.5.2   
#  34 hwriter          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.3.2   
#  35 ini              /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.1   
#  36 IRanges          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.16.0  
#  37 jsonlite         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.6     
#  38 labeling         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3     
#  39 lambda.r         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.2.4   
#  40 latticeExtra     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.6-28  
#  41 lazyeval         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.2.2   
#  42 lifecycle        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.1.0   
#  43 magrittr         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.5     
#  44 matrixStats      /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.55.0  
#  45 munsell          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.5.0   
#  46 packrat          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.5.0   
#  47 pillar           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4.2   
#  48 pkgconfig        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.0.3   
#  49 plogr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.2.0   
#  50 plyr             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.8.5   
#  51 purrr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.3   
#  52 R6               /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.4.1   
#  53 RColorBrewer     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.1-2   
#  54 Rcpp             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.0.3   
#  55 RcppParallel     /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 4.4.4   
#  56 RcppTOML         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.1.6   
#  57 RCurl            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.95-4.…
#  58 readr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.3.1   
#  59 reshape2         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4.3   
#  60 rlang            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.4.2   
#  61 Rsamtools        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.34.1  
#  62 S4Vectors        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.20.1  
#  63 scales           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.1.0   
#  64 ShortRead        /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.40.0  
#  65 snow             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.4-3   
#  66 stringi          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4.3   
#  67 stringr          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.4.0   
#  68 SummarizedExper… /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.12.0  
#  69 tibble           /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.1.3   
#  70 tidyselect       /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.2.5   
#  71 utf8             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.1.4   
#  72 vctrs            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.2.0   
#  73 viridisLite      /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.3.0   
#  74 withr            /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.1.2   
#  75 XVector          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.22.0  
#  76 yaml             /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 2.2.0   
#  77 zeallot          /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 0.1.0   
#  78 zlibbioc         /home/btupper/edna/edna-dada2/packrat/lib/x86_64-… 1.28.0  
#  79 base             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  80 class            /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 7.3-15  
#  81 codetools        /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 0.2-16  
#  82 compiler         /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  83 datasets         /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  84 foreign          /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 0.8-71  
#  85 graphics         /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  86 grDevices        /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  87 grid             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  88 lattice          /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 0.20-38 
#  89 MASS             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 7.3-51.3
#  90 Matrix           /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 1.2-17  
#  91 methods          /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  92 mgcv             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 1.8-28  
#  93 nlme             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.1-137 
#  94 nnet             /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 7.3-12  
#  95 parallel         /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  96 splines          /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  97 stats            /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  98 stats4           /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
#  99 survival         /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 2.44-1.1
# 100 tcltk            /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
# 101 tools            /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1   
# 102 utils            /home/btupper/edna/edna-dada2/packrat/lib-R/x86_6… 3.5.1  
```