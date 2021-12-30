Projet CC2 paillomavirus
================

Margaux BONNARDOT M1 MFA

# Problématique :

Lors d’analyses de microbiote, l’utilisation d’Operational Taxonomic
Units(OTU) qui consite à un clustering des sequences qui ont plus de 97%
de similarité a longtemps été la methode de regourpement la plus
couragment utilisée. On utilise ce clustering en metabarrecoding pour
regrouper les individus phylogénétiquement proches au seins du
microbiote. Mais la notion d’espèce est difficile à qualifier en terme
de pourcentage de similitude pour une certaine sequence. D’autres
methodes on ainsi été developpées. On peut notamment citer DADA2 qui ne
repose pas sur le regroupement par cluster basé sur la similarité. Cette
méthode utilise des Amplicon sequencing variants (ASV). Cette méthode
est basée sur une approche probabiliste et permet la détection et la
correction des erreurs de séquençage et dépendant du jeu de données.
Cela permet contrairement aux méthodes de novo et closed-ref de
connaître l’abondance réelle car les erreurs de séquençage sont
éliminées et on se retrouve avec seulement les séquences différentes
d’un point de vue biologique.

Nous souhaitons regarder l’impact au niveau de l’alpha et la beta
diversité suivant la méthode utilisé. Pour cela nous utilisont les
données de l’article “Depiction of Vaginal Microbiota in Women With
High-Risk Human Papillomavirus Infection” qui utilise le clustering à
97% et regarder si une difference significative peut être observable au
niveau des indices de Shannon, de Chao1 et d’Observed species ainsi que
de principal coordinate. Les données que nous analysons ici sont des
séquences d’amplicons Illumina Miseq 2x250 et il s’agit du séquençage de
la zone hypervariable V3-V4 de l’ARNr 16S.

# Traitement des données :

## Packages utilisées :

``` r
library("knitr")
library("BiocStyle")
```

``` r
.cran_packages <- c("ggplot2", "gridExtra", "devtools")

.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ## Loading required package: ggplot2

    ## Loading required package: gridExtra

    ## Loading required package: devtools

    ## Loading required package: usethis

    ## Loading required package: dada2

    ## Loading required package: Rcpp

    ## Loading required package: phyloseq

    ## Loading required package: DECIPHER

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

    ## Loading required package: phangorn

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

    ##   ggplot2 gridExtra  devtools     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools", "reshape2", "PMA", "structSSI", "ade4", "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")

# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'structSSI' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

``` r
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
library(ggrepel)
library(phyloseq)
library(ggplot2)
library(dada2)
library(DECIPHER)
library(phangorn)
```

## Traitement des données brutes

``` r
set.seed(120)
miseq_path <- "/home/rstudio/HPVsegfastq" 
list.files(miseq_path)
```

    ##   [1] "filtered"              "SRR9611359_1.fastq.gz" "SRR9611359_2.fastq.gz"
    ##   [4] "SRR9611360_1.fastq.gz" "SRR9611360_2.fastq.gz" "SRR9611361_1.fastq.gz"
    ##   [7] "SRR9611361_2.fastq.gz" "SRR9611362_1.fastq.gz" "SRR9611362_2.fastq.gz"
    ##  [10] "SRR9611363_1.fastq.gz" "SRR9611363_2.fastq.gz" "SRR9611364_1.fastq.gz"
    ##  [13] "SRR9611364_2.fastq.gz" "SRR9611365_1.fastq.gz" "SRR9611365_2.fastq.gz"
    ##  [16] "SRR9611366_1.fastq.gz" "SRR9611366_2.fastq.gz" "SRR9611367_1.fastq.gz"
    ##  [19] "SRR9611367_2.fastq.gz" "SRR9611368_1.fastq.gz" "SRR9611368_2.fastq.gz"
    ##  [22] "SRR9611369_1.fastq.gz" "SRR9611369_2.fastq.gz" "SRR9611370_1.fastq.gz"
    ##  [25] "SRR9611370_2.fastq.gz" "SRR9611371_1.fastq.gz" "SRR9611371_2.fastq.gz"
    ##  [28] "SRR9611372_1.fastq.gz" "SRR9611372_2.fastq.gz" "SRR9611373_1.fastq.gz"
    ##  [31] "SRR9611373_2.fastq.gz" "SRR9611374_1.fastq.gz" "SRR9611374_2.fastq.gz"
    ##  [34] "SRR9611375_1.fastq.gz" "SRR9611375_2.fastq.gz" "SRR9611376_1.fastq.gz"
    ##  [37] "SRR9611376_2.fastq.gz" "SRR9611377_1.fastq.gz" "SRR9611377_2.fastq.gz"
    ##  [40] "SRR9611378_1.fastq.gz" "SRR9611378_2.fastq.gz" "SRR9611379_1.fastq.gz"
    ##  [43] "SRR9611379_2.fastq.gz" "SRR9611380_1.fastq.gz" "SRR9611380_2.fastq.gz"
    ##  [46] "SRR9611381_1.fastq.gz" "SRR9611381_2.fastq.gz" "SRR9611382_1.fastq.gz"
    ##  [49] "SRR9611382_2.fastq.gz" "SRR9611383_1.fastq.gz" "SRR9611383_2.fastq.gz"
    ##  [52] "SRR9611384_1.fastq.gz" "SRR9611384_2.fastq.gz" "SRR9611385_1.fastq.gz"
    ##  [55] "SRR9611385_2.fastq.gz" "SRR9611386_1.fastq.gz" "SRR9611386_2.fastq.gz"
    ##  [58] "SRR9611387_1.fastq.gz" "SRR9611387_2.fastq.gz" "SRR9611388_1.fastq.gz"
    ##  [61] "SRR9611388_2.fastq.gz" "SRR9611389_1.fastq.gz" "SRR9611389_2.fastq.gz"
    ##  [64] "SRR9611390_1.fastq.gz" "SRR9611390_2.fastq.gz" "SRR9611391_1.fastq.gz"
    ##  [67] "SRR9611391_2.fastq.gz" "SRR9611392_1.fastq.gz" "SRR9611392_2.fastq.gz"
    ##  [70] "SRR9611393_1.fastq.gz" "SRR9611393_2.fastq.gz" "SRR9611394_1.fastq.gz"
    ##  [73] "SRR9611394_2.fastq.gz" "SRR9611395_1.fastq.gz" "SRR9611395_2.fastq.gz"
    ##  [76] "SRR9611396_1.fastq.gz" "SRR9611396_2.fastq.gz" "SRR9611397_1.fastq.gz"
    ##  [79] "SRR9611397_2.fastq.gz" "SRR9611398_1.fastq.gz" "SRR9611398_2.fastq.gz"
    ##  [82] "SRR9611399_1.fastq.gz" "SRR9611399_2.fastq.gz" "SRR9611400_1.fastq.gz"
    ##  [85] "SRR9611400_2.fastq.gz" "SRR9611401_1.fastq.gz" "SRR9611401_2.fastq.gz"
    ##  [88] "SRR9611402_1.fastq.gz" "SRR9611402_2.fastq.gz" "SRR9611403_1.fastq.gz"
    ##  [91] "SRR9611403_2.fastq.gz" "SRR9611404_1.fastq.gz" "SRR9611404_2.fastq.gz"
    ##  [94] "SRR9611405_1.fastq.gz" "SRR9611405_2.fastq.gz" "SRR9611406_1.fastq.gz"
    ##  [97] "SRR9611406_2.fastq.gz" "SRR9611407_1.fastq.gz" "SRR9611407_2.fastq.gz"
    ## [100] "SRR9611408_1.fastq.gz" "SRR9611408_2.fastq.gz" "SRR9611409_1.fastq.gz"
    ## [103] "SRR9611409_2.fastq.gz" "SRR9611410_1.fastq.gz" "SRR9611410_2.fastq.gz"
    ## [106] "SRR9611411_1.fastq.gz" "SRR9611411_2.fastq.gz" "SRR9611412_1.fastq.gz"
    ## [109] "SRR9611412_2.fastq.gz" "SRR9611413_1.fastq.gz" "SRR9611413_2.fastq.gz"
    ## [112] "SRR9611414_1.fastq.gz" "SRR9611414_2.fastq.gz" "SRR9611415_1.fastq.gz"
    ## [115] "SRR9611415_2.fastq.gz" "SRR9611416_1.fastq.gz" "SRR9611416_2.fastq.gz"
    ## [118] "SRR9611417_1.fastq.gz" "SRR9611417_2.fastq.gz" "SRR9611418_1.fastq.gz"
    ## [121] "SRR9611418_2.fastq.gz"

### Filter and Trim

Je creer 2 objets qui vont contenir les données des fichiers fastq, les
read 1 sont pour les forward et les read 2 pour les reverse. La fonction
sort permet de classer par ordre alphabetique. On specifit le path vers
fnFs et fnRs avec file.path. Pour verifier que la commande precedante à
bien marché, on affiches les 3 premiers éléments de fnFsHPV et fnRsHPV.
On obtient:

``` r
fnFs <- sort(list.files(miseq_path, pattern="_1.fastq.gz")) 
fnRs <- sort(list.files(miseq_path, pattern="_2.fastq.gz")) 

sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)


fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

fnFs[1:3]
```

    ## [1] "/home/rstudio/HPVsegfastq/SRR9611359_1.fastq.gz"
    ## [2] "/home/rstudio/HPVsegfastq/SRR9611360_1.fastq.gz"
    ## [3] "/home/rstudio/HPVsegfastq/SRR9611361_1.fastq.gz"

``` r
fnRs[1:3]
```

    ## [1] "/home/rstudio/HPVsegfastq/SRR9611359_2.fastq.gz"
    ## [2] "/home/rstudio/HPVsegfastq/SRR9611360_2.fastq.gz"
    ## [3] "/home/rstudio/HPVsegfastq/SRR9611361_2.fastq.gz"

``` r
# Pour les fowards:
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Pour les reverse:
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

Cette etape permet de regarder pour les 2 premiers echantillions la
proportion des reads qui sont à un certain quality score pour une base
donnée. On regarde donc ce qu’il faut trim, on considère qu’il faut
avoir des reads avec un quality score de au moin 30.

Pour trim mes sequences je trim la partie de droite à la longueur de 240
bp

``` r
filt_path <- file.path(miseq_path, "filtered") 

if(!file_test("-d", filt_path)) 
  dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fasta.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2),truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)

head(out)
```

    ##                       reads.in reads.out
    ## SRR9611359_1.fastq.gz    40550     36736
    ## SRR9611360_1.fastq.gz    43567     39713
    ## SRR9611361_1.fastq.gz    33032     28214
    ## SRR9611362_1.fastq.gz    54486     43700
    ## SRR9611363_1.fastq.gz    44036     35129
    ## SRR9611364_1.fastq.gz    55035     46545

#### Dereplication

``` r
  derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611359_F_filt.fasta.gz

    ## Encountered 6543 unique sequences from 36736 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611360_F_filt.fasta.gz

    ## Encountered 7373 unique sequences from 39713 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611361_F_filt.fasta.gz

    ## Encountered 8522 unique sequences from 28214 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611362_F_filt.fasta.gz

    ## Encountered 20644 unique sequences from 43700 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611363_F_filt.fasta.gz

    ## Encountered 13296 unique sequences from 35129 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611364_F_filt.fasta.gz

    ## Encountered 18883 unique sequences from 46545 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611365_F_filt.fasta.gz

    ## Encountered 16780 unique sequences from 41125 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611366_F_filt.fasta.gz

    ## Encountered 8822 unique sequences from 29068 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611367_F_filt.fasta.gz

    ## Encountered 10138 unique sequences from 48836 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611368_F_filt.fasta.gz

    ## Encountered 7252 unique sequences from 38695 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611369_F_filt.fasta.gz

    ## Encountered 6929 unique sequences from 40412 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611370_F_filt.fasta.gz

    ## Encountered 9518 unique sequences from 55227 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611371_F_filt.fasta.gz

    ## Encountered 7798 unique sequences from 42638 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611372_F_filt.fasta.gz

    ## Encountered 8997 unique sequences from 45520 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611373_F_filt.fasta.gz

    ## Encountered 11640 unique sequences from 60508 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611374_F_filt.fasta.gz

    ## Encountered 8319 unique sequences from 52186 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611375_F_filt.fasta.gz

    ## Encountered 9634 unique sequences from 51031 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611376_F_filt.fasta.gz

    ## Encountered 9169 unique sequences from 48157 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611377_F_filt.fasta.gz

    ## Encountered 6859 unique sequences from 37435 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611378_F_filt.fasta.gz

    ## Encountered 8019 unique sequences from 46201 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611379_F_filt.fasta.gz

    ## Encountered 9785 unique sequences from 50515 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611380_F_filt.fasta.gz

    ## Encountered 8797 unique sequences from 49834 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611381_F_filt.fasta.gz

    ## Encountered 8002 unique sequences from 44996 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611382_F_filt.fasta.gz

    ## Encountered 8245 unique sequences from 46762 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611383_F_filt.fasta.gz

    ## Encountered 14030 unique sequences from 54583 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611384_F_filt.fasta.gz

    ## Encountered 8599 unique sequences from 51905 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611385_F_filt.fasta.gz

    ## Encountered 11292 unique sequences from 48684 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611386_F_filt.fasta.gz

    ## Encountered 11835 unique sequences from 61364 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611387_F_filt.fasta.gz

    ## Encountered 9202 unique sequences from 50491 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611388_F_filt.fasta.gz

    ## Encountered 8837 unique sequences from 48133 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611389_F_filt.fasta.gz

    ## Encountered 5922 unique sequences from 27221 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611390_F_filt.fasta.gz

    ## Encountered 9375 unique sequences from 29482 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611391_F_filt.fasta.gz

    ## Encountered 20420 unique sequences from 48275 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611392_F_filt.fasta.gz

    ## Encountered 11243 unique sequences from 54841 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611393_F_filt.fasta.gz

    ## Encountered 8468 unique sequences from 40826 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611394_F_filt.fasta.gz

    ## Encountered 7701 unique sequences from 47260 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611395_F_filt.fasta.gz

    ## Encountered 11127 unique sequences from 57153 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611396_F_filt.fasta.gz

    ## Encountered 9544 unique sequences from 56286 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611397_F_filt.fasta.gz

    ## Encountered 11238 unique sequences from 54383 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611398_F_filt.fasta.gz

    ## Encountered 7797 unique sequences from 45525 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611399_F_filt.fasta.gz

    ## Encountered 9556 unique sequences from 49024 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611400_F_filt.fasta.gz

    ## Encountered 18102 unique sequences from 89619 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611401_F_filt.fasta.gz

    ## Encountered 10564 unique sequences from 54726 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611402_F_filt.fasta.gz

    ## Encountered 9202 unique sequences from 52970 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611403_F_filt.fasta.gz

    ## Encountered 7561 unique sequences from 39241 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611404_F_filt.fasta.gz

    ## Encountered 8363 unique sequences from 42587 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611405_F_filt.fasta.gz

    ## Encountered 6666 unique sequences from 41072 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611406_F_filt.fasta.gz

    ## Encountered 7290 unique sequences from 39099 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611407_F_filt.fasta.gz

    ## Encountered 7075 unique sequences from 43134 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611408_F_filt.fasta.gz

    ## Encountered 7705 unique sequences from 43207 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611409_F_filt.fasta.gz

    ## Encountered 7256 unique sequences from 42050 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611410_F_filt.fasta.gz

    ## Encountered 15737 unique sequences from 106910 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611411_F_filt.fasta.gz

    ## Encountered 10939 unique sequences from 58669 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611412_F_filt.fasta.gz

    ## Encountered 8916 unique sequences from 55227 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611413_F_filt.fasta.gz

    ## Encountered 8304 unique sequences from 45988 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611414_F_filt.fasta.gz

    ## Encountered 7043 unique sequences from 38089 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611415_F_filt.fasta.gz

    ## Encountered 12051 unique sequences from 51054 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611416_F_filt.fasta.gz

    ## Encountered 8073 unique sequences from 44568 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611417_F_filt.fasta.gz

    ## Encountered 12841 unique sequences from 72580 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611418_F_filt.fasta.gz

    ## Encountered 17503 unique sequences from 100127 total sequences read.

``` r
  derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611359_R_filt.fastq.gz

    ## Encountered 7096 unique sequences from 36736 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611360_R_filt.fastq.gz

    ## Encountered 8092 unique sequences from 39713 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611361_R_filt.fastq.gz

    ## Encountered 8540 unique sequences from 28214 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611362_R_filt.fastq.gz

    ## Encountered 12187 unique sequences from 43700 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611363_R_filt.fastq.gz

    ## Encountered 12498 unique sequences from 35129 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611364_R_filt.fastq.gz

    ## Encountered 12954 unique sequences from 46545 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611365_R_filt.fastq.gz

    ## Encountered 13046 unique sequences from 41125 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611366_R_filt.fastq.gz

    ## Encountered 9565 unique sequences from 29068 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611367_R_filt.fastq.gz

    ## Encountered 10658 unique sequences from 48836 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611368_R_filt.fastq.gz

    ## Encountered 8326 unique sequences from 38695 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611369_R_filt.fastq.gz

    ## Encountered 7820 unique sequences from 40412 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611370_R_filt.fastq.gz

    ## Encountered 9620 unique sequences from 55227 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611371_R_filt.fastq.gz

    ## Encountered 7851 unique sequences from 42638 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611372_R_filt.fastq.gz

    ## Encountered 10962 unique sequences from 45520 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611373_R_filt.fastq.gz

    ## Encountered 11526 unique sequences from 60508 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611374_R_filt.fastq.gz

    ## Encountered 8447 unique sequences from 52186 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611375_R_filt.fastq.gz

    ## Encountered 10561 unique sequences from 51031 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611376_R_filt.fastq.gz

    ## Encountered 9536 unique sequences from 48157 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611377_R_filt.fastq.gz

    ## Encountered 7835 unique sequences from 37435 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611378_R_filt.fastq.gz

    ## Encountered 8044 unique sequences from 46201 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611379_R_filt.fastq.gz

    ## Encountered 10909 unique sequences from 50515 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611380_R_filt.fastq.gz

    ## Encountered 9315 unique sequences from 49834 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611381_R_filt.fastq.gz

    ## Encountered 8704 unique sequences from 44996 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611382_R_filt.fastq.gz

    ## Encountered 9647 unique sequences from 46762 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611383_R_filt.fastq.gz

    ## Encountered 15522 unique sequences from 54583 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611384_R_filt.fastq.gz

    ## Encountered 9464 unique sequences from 51905 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611385_R_filt.fastq.gz

    ## Encountered 12105 unique sequences from 48684 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611386_R_filt.fastq.gz

    ## Encountered 13288 unique sequences from 61364 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611387_R_filt.fastq.gz

    ## Encountered 9230 unique sequences from 50491 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611388_R_filt.fastq.gz

    ## Encountered 9702 unique sequences from 48133 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611389_R_filt.fastq.gz

    ## Encountered 5924 unique sequences from 27221 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611390_R_filt.fastq.gz

    ## Encountered 9345 unique sequences from 29482 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611391_R_filt.fastq.gz

    ## Encountered 12794 unique sequences from 48275 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611392_R_filt.fastq.gz

    ## Encountered 14137 unique sequences from 54841 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611393_R_filt.fastq.gz

    ## Encountered 9755 unique sequences from 40826 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611394_R_filt.fastq.gz

    ## Encountered 8879 unique sequences from 47260 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611395_R_filt.fastq.gz

    ## Encountered 5902 unique sequences from 57153 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611396_R_filt.fastq.gz

    ## Encountered 9675 unique sequences from 56286 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611397_R_filt.fastq.gz

    ## Encountered 6114 unique sequences from 54383 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611398_R_filt.fastq.gz

    ## Encountered 10335 unique sequences from 45525 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611399_R_filt.fastq.gz

    ## Encountered 9714 unique sequences from 49024 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611400_R_filt.fastq.gz

    ## Encountered 16243 unique sequences from 89619 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611401_R_filt.fastq.gz

    ## Encountered 11852 unique sequences from 54726 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611402_R_filt.fastq.gz

    ## Encountered 10723 unique sequences from 52970 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611403_R_filt.fastq.gz

    ## Encountered 10208 unique sequences from 39241 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611404_R_filt.fastq.gz

    ## Encountered 9998 unique sequences from 42587 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611405_R_filt.fastq.gz

    ## Encountered 8630 unique sequences from 41072 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611406_R_filt.fastq.gz

    ## Encountered 7997 unique sequences from 39099 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611407_R_filt.fastq.gz

    ## Encountered 8257 unique sequences from 43134 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611408_R_filt.fastq.gz

    ## Encountered 9173 unique sequences from 43207 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611409_R_filt.fastq.gz

    ## Encountered 7518 unique sequences from 42050 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611410_R_filt.fastq.gz

    ## Encountered 14066 unique sequences from 106910 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611411_R_filt.fastq.gz

    ## Encountered 12740 unique sequences from 58669 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611412_R_filt.fastq.gz

    ## Encountered 8991 unique sequences from 55227 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611413_R_filt.fastq.gz

    ## Encountered 9073 unique sequences from 45988 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611414_R_filt.fastq.gz

    ## Encountered 8218 unique sequences from 38089 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611415_R_filt.fastq.gz

    ## Encountered 14056 unique sequences from 51054 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611416_R_filt.fastq.gz

    ## Encountered 9695 unique sequences from 44568 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611417_R_filt.fastq.gz

    ## Encountered 15299 unique sequences from 72580 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/HPVsegfastq/filtered/SRR9611418_R_filt.fastq.gz

    ## Encountered 18744 unique sequences from 100127 total sequences read.

``` r
    names(derepFs) <- sampleNames
    names(derepRs) <- sampleNames
```

``` r
  errF <- learnErrors(filtFs, multithread=TRUE) 
```

    ## 102761520 total bases in 428173 reads from 11 samples will be used for learning the error rates.

``` r
  errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 101130560 total bases in 632066 reads from 15 samples will be used for learning the error rates.

``` r
  plotErrors(errF)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
  plotErrors(errR)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

Pooling improves the detection of rare variants

``` r
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 36736 reads in 6543 unique sequences.
    ## Sample 2 - 39713 reads in 7373 unique sequences.
    ## Sample 3 - 28214 reads in 8522 unique sequences.
    ## Sample 4 - 43700 reads in 20644 unique sequences.
    ## Sample 5 - 35129 reads in 13296 unique sequences.
    ## Sample 6 - 46545 reads in 18883 unique sequences.
    ## Sample 7 - 41125 reads in 16780 unique sequences.
    ## Sample 8 - 29068 reads in 8822 unique sequences.
    ## Sample 9 - 48836 reads in 10138 unique sequences.
    ## Sample 10 - 38695 reads in 7252 unique sequences.
    ## Sample 11 - 40412 reads in 6929 unique sequences.
    ## Sample 12 - 55227 reads in 9518 unique sequences.
    ## Sample 13 - 42638 reads in 7798 unique sequences.
    ## Sample 14 - 45520 reads in 8997 unique sequences.
    ## Sample 15 - 60508 reads in 11640 unique sequences.
    ## Sample 16 - 52186 reads in 8319 unique sequences.
    ## Sample 17 - 51031 reads in 9634 unique sequences.
    ## Sample 18 - 48157 reads in 9169 unique sequences.
    ## Sample 19 - 37435 reads in 6859 unique sequences.
    ## Sample 20 - 46201 reads in 8019 unique sequences.
    ## Sample 21 - 50515 reads in 9785 unique sequences.
    ## Sample 22 - 49834 reads in 8797 unique sequences.
    ## Sample 23 - 44996 reads in 8002 unique sequences.
    ## Sample 24 - 46762 reads in 8245 unique sequences.
    ## Sample 25 - 54583 reads in 14030 unique sequences.
    ## Sample 26 - 51905 reads in 8599 unique sequences.
    ## Sample 27 - 48684 reads in 11292 unique sequences.
    ## Sample 28 - 61364 reads in 11835 unique sequences.
    ## Sample 29 - 50491 reads in 9202 unique sequences.
    ## Sample 30 - 48133 reads in 8837 unique sequences.
    ## Sample 31 - 27221 reads in 5922 unique sequences.
    ## Sample 32 - 29482 reads in 9375 unique sequences.
    ## Sample 33 - 48275 reads in 20420 unique sequences.
    ## Sample 34 - 54841 reads in 11243 unique sequences.
    ## Sample 35 - 40826 reads in 8468 unique sequences.
    ## Sample 36 - 47260 reads in 7701 unique sequences.
    ## Sample 37 - 57153 reads in 11127 unique sequences.
    ## Sample 38 - 56286 reads in 9544 unique sequences.
    ## Sample 39 - 54383 reads in 11238 unique sequences.
    ## Sample 40 - 45525 reads in 7797 unique sequences.
    ## Sample 41 - 49024 reads in 9556 unique sequences.
    ## Sample 42 - 89619 reads in 18102 unique sequences.
    ## Sample 43 - 54726 reads in 10564 unique sequences.
    ## Sample 44 - 52970 reads in 9202 unique sequences.
    ## Sample 45 - 39241 reads in 7561 unique sequences.
    ## Sample 46 - 42587 reads in 8363 unique sequences.
    ## Sample 47 - 41072 reads in 6666 unique sequences.
    ## Sample 48 - 39099 reads in 7290 unique sequences.
    ## Sample 49 - 43134 reads in 7075 unique sequences.
    ## Sample 50 - 43207 reads in 7705 unique sequences.
    ## Sample 51 - 42050 reads in 7256 unique sequences.
    ## Sample 52 - 106910 reads in 15737 unique sequences.
    ## Sample 53 - 58669 reads in 10939 unique sequences.
    ## Sample 54 - 55227 reads in 8916 unique sequences.
    ## Sample 55 - 45988 reads in 8304 unique sequences.
    ## Sample 56 - 38089 reads in 7043 unique sequences.
    ## Sample 57 - 51054 reads in 12051 unique sequences.
    ## Sample 58 - 44568 reads in 8073 unique sequences.
    ## Sample 59 - 72580 reads in 12841 unique sequences.
    ## Sample 60 - 100127 reads in 17503 unique sequences.

``` r
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 36736 reads in 7096 unique sequences.
    ## Sample 2 - 39713 reads in 8092 unique sequences.
    ## Sample 3 - 28214 reads in 8540 unique sequences.
    ## Sample 4 - 43700 reads in 12187 unique sequences.
    ## Sample 5 - 35129 reads in 12498 unique sequences.
    ## Sample 6 - 46545 reads in 12954 unique sequences.
    ## Sample 7 - 41125 reads in 13046 unique sequences.
    ## Sample 8 - 29068 reads in 9565 unique sequences.
    ## Sample 9 - 48836 reads in 10658 unique sequences.
    ## Sample 10 - 38695 reads in 8326 unique sequences.
    ## Sample 11 - 40412 reads in 7820 unique sequences.
    ## Sample 12 - 55227 reads in 9620 unique sequences.
    ## Sample 13 - 42638 reads in 7851 unique sequences.
    ## Sample 14 - 45520 reads in 10962 unique sequences.
    ## Sample 15 - 60508 reads in 11526 unique sequences.
    ## Sample 16 - 52186 reads in 8447 unique sequences.
    ## Sample 17 - 51031 reads in 10561 unique sequences.
    ## Sample 18 - 48157 reads in 9536 unique sequences.
    ## Sample 19 - 37435 reads in 7835 unique sequences.
    ## Sample 20 - 46201 reads in 8044 unique sequences.
    ## Sample 21 - 50515 reads in 10909 unique sequences.
    ## Sample 22 - 49834 reads in 9315 unique sequences.
    ## Sample 23 - 44996 reads in 8704 unique sequences.
    ## Sample 24 - 46762 reads in 9647 unique sequences.
    ## Sample 25 - 54583 reads in 15522 unique sequences.
    ## Sample 26 - 51905 reads in 9464 unique sequences.
    ## Sample 27 - 48684 reads in 12105 unique sequences.
    ## Sample 28 - 61364 reads in 13288 unique sequences.
    ## Sample 29 - 50491 reads in 9230 unique sequences.
    ## Sample 30 - 48133 reads in 9702 unique sequences.
    ## Sample 31 - 27221 reads in 5924 unique sequences.
    ## Sample 32 - 29482 reads in 9345 unique sequences.
    ## Sample 33 - 48275 reads in 12794 unique sequences.
    ## Sample 34 - 54841 reads in 14137 unique sequences.
    ## Sample 35 - 40826 reads in 9755 unique sequences.
    ## Sample 36 - 47260 reads in 8879 unique sequences.
    ## Sample 37 - 57153 reads in 5902 unique sequences.
    ## Sample 38 - 56286 reads in 9675 unique sequences.
    ## Sample 39 - 54383 reads in 6114 unique sequences.
    ## Sample 40 - 45525 reads in 10335 unique sequences.
    ## Sample 41 - 49024 reads in 9714 unique sequences.
    ## Sample 42 - 89619 reads in 16243 unique sequences.
    ## Sample 43 - 54726 reads in 11852 unique sequences.
    ## Sample 44 - 52970 reads in 10723 unique sequences.
    ## Sample 45 - 39241 reads in 10208 unique sequences.
    ## Sample 46 - 42587 reads in 9998 unique sequences.
    ## Sample 47 - 41072 reads in 8630 unique sequences.
    ## Sample 48 - 39099 reads in 7997 unique sequences.
    ## Sample 49 - 43134 reads in 8257 unique sequences.
    ## Sample 50 - 43207 reads in 9173 unique sequences.
    ## Sample 51 - 42050 reads in 7518 unique sequences.
    ## Sample 52 - 106910 reads in 14066 unique sequences.
    ## Sample 53 - 58669 reads in 12740 unique sequences.
    ## Sample 54 - 55227 reads in 8991 unique sequences.
    ## Sample 55 - 45988 reads in 9073 unique sequences.
    ## Sample 56 - 38089 reads in 8218 unique sequences.
    ## Sample 57 - 51054 reads in 14056 unique sequences.
    ## Sample 58 - 44568 reads in 9695 unique sequences.
    ## Sample 59 - 72580 reads in 15299 unique sequences.
    ## Sample 60 - 100127 reads in 18744 unique sequences.

Inspecting the dada-class object returned by dada:

``` r
  dadaFs[[1]]  
```

    ## dada-class: object describing DADA2 denoising results
    ## 5 sequence variants were inferred from 6543 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
  dadaRs[[1]]  
```

    ## dada-class: object describing DADA2 denoising results
    ## 21 sequence variants were inferred from 7096 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

### Construct sequence table and remove chimeras

``` r
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

  seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
  table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 295 304 305 313 
    ##   1   5   1   6

``` r
  seqtabNoC <- removeBimeraDenovo(seqtabAll)
  seqtabNoC
```

    ##            ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                 6
    ## SRR9611365                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                 0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                 0
    ##            ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT
    ## SRR9611359                                                                                                                                                                                                                                                                                                       0
    ## SRR9611360                                                                                                                                                                                                                                                                                                       0
    ## SRR9611361                                                                                                                                                                                                                                                                                                       0
    ## SRR9611362                                                                                                                                                                                                                                                                                                       0
    ## SRR9611363                                                                                                                                                                                                                                                                                                       0
    ## SRR9611364                                                                                                                                                                                                                                                                                                       0
    ## SRR9611365                                                                                                                                                                                                                                                                                                       0
    ## SRR9611366                                                                                                                                                                                                                                                                                                       5
    ## SRR9611367                                                                                                                                                                                                                                                                                                       0
    ## SRR9611368                                                                                                                                                                                                                                                                                                       0
    ## SRR9611369                                                                                                                                                                                                                                                                                                       0
    ## SRR9611370                                                                                                                                                                                                                                                                                                       0
    ## SRR9611371                                                                                                                                                                                                                                                                                                       0
    ## SRR9611372                                                                                                                                                                                                                                                                                                       0
    ## SRR9611373                                                                                                                                                                                                                                                                                                       0
    ## SRR9611374                                                                                                                                                                                                                                                                                                       0
    ## SRR9611375                                                                                                                                                                                                                                                                                                       0
    ## SRR9611376                                                                                                                                                                                                                                                                                                       0
    ## SRR9611377                                                                                                                                                                                                                                                                                                       0
    ## SRR9611378                                                                                                                                                                                                                                                                                                       0
    ## SRR9611379                                                                                                                                                                                                                                                                                                       0
    ## SRR9611380                                                                                                                                                                                                                                                                                                       0
    ## SRR9611381                                                                                                                                                                                                                                                                                                       0
    ## SRR9611382                                                                                                                                                                                                                                                                                                       0
    ## SRR9611383                                                                                                                                                                                                                                                                                                       0
    ## SRR9611384                                                                                                                                                                                                                                                                                                       0
    ## SRR9611385                                                                                                                                                                                                                                                                                                       0
    ## SRR9611386                                                                                                                                                                                                                                                                                                       0
    ## SRR9611387                                                                                                                                                                                                                                                                                                       0
    ## SRR9611388                                                                                                                                                                                                                                                                                                       0
    ## SRR9611389                                                                                                                                                                                                                                                                                                       0
    ## SRR9611390                                                                                                                                                                                                                                                                                                       0
    ## SRR9611391                                                                                                                                                                                                                                                                                                       0
    ## SRR9611392                                                                                                                                                                                                                                                                                                       0
    ## SRR9611393                                                                                                                                                                                                                                                                                                       0
    ## SRR9611394                                                                                                                                                                                                                                                                                                       0
    ## SRR9611395                                                                                                                                                                                                                                                                                                       0
    ## SRR9611396                                                                                                                                                                                                                                                                                                       0
    ## SRR9611397                                                                                                                                                                                                                                                                                                       0
    ## SRR9611398                                                                                                                                                                                                                                                                                                       0
    ## SRR9611399                                                                                                                                                                                                                                                                                                       0
    ## SRR9611400                                                                                                                                                                                                                                                                                                       0
    ## SRR9611401                                                                                                                                                                                                                                                                                                       0
    ## SRR9611402                                                                                                                                                                                                                                                                                                       0
    ## SRR9611403                                                                                                                                                                                                                                                                                                       0
    ## SRR9611404                                                                                                                                                                                                                                                                                                       0
    ## SRR9611405                                                                                                                                                                                                                                                                                                       0
    ## SRR9611406                                                                                                                                                                                                                                                                                                       0
    ## SRR9611407                                                                                                                                                                                                                                                                                                       0
    ## SRR9611408                                                                                                                                                                                                                                                                                                       0
    ## SRR9611409                                                                                                                                                                                                                                                                                                       0
    ## SRR9611410                                                                                                                                                                                                                                                                                                       0
    ## SRR9611411                                                                                                                                                                                                                                                                                                       0
    ## SRR9611412                                                                                                                                                                                                                                                                                                       0
    ## SRR9611413                                                                                                                                                                                                                                                                                                       0
    ## SRR9611414                                                                                                                                                                                                                                                                                                       0
    ## SRR9611415                                                                                                                                                                                                                                                                                                       0
    ## SRR9611416                                                                                                                                                                                                                                                                                                       0
    ## SRR9611417                                                                                                                                                                                                                                                                                                       0
    ## SRR9611418                                                                                                                                                                                                                                                                                                       0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         2
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0
    ##            TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT
    ## SRR9611359                                                                                                                                                                                                                                                                                                                0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                0
    ## SRR9611365                                                                                                                                                                                                                                                                                                                0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                2
    ## SRR9611401                                                                                                                                                                                                                                                                                                                0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                0
    ##            TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT
    ## SRR9611359                                                                                                                                                                                                                                                                                                                0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                0
    ## SRR9611365                                                                                                                                                                                                                                                                                                                0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                2
    ## SRR9611401                                                                                                                                                                                                                                                                                                                0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                0
    ##            TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT
    ## SRR9611359                                                                                                                                                                                                                                                                                                                0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                0
    ## SRR9611365                                                                                                                                                                                                                                                                                                                0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                2
    ## SRR9611401                                                                                                                                                                                                                                                                                                                0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         1
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         1
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         1
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         1
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0
    ##            ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG
    ## SRR9611359                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611360                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611361                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611362                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611363                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611364                                                                                                                                                                                                                                                                                                                         1
    ## SRR9611365                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611366                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611367                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611368                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611369                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611370                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611371                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611372                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611373                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611374                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611375                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611376                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611377                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611378                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611379                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611380                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611381                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611382                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611383                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611384                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611385                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611386                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611387                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611388                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611389                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611390                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611391                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611392                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611393                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611394                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611395                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611396                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611397                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611398                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611399                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611400                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611401                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611402                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611403                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611404                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611405                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611406                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611407                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611408                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611409                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611410                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611411                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611412                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611413                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611414                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611415                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611416                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611417                                                                                                                                                                                                                                                                                                                         0
    ## SRR9611418                                                                                                                                                                                                                                                                                                                         0

### Assign taxonomy

``` r
  fastaRef <- "/home/rstudio/silva_nr99_v138.1_train_set.fa"
  taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE) 
  unname(head(taxTab))
```

    ##      [,1]       [,2]             [,3]                  [,4]             
    ## [1,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rickettsiales"  
    ## [2,] "Bacteria" "Proteobacteria" "Alphaproteobacteria" "Rickettsiales"  
    ## [3,] "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"
    ## [4,] "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"
    ## [5,] "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"
    ## [6,] "Bacteria" "Firmicutes"     "Bacilli"             "Lactobacillales"
    ##      [,5]               [,6]           
    ## [1,] "Mitochondria"     NA             
    ## [2,] "Mitochondria"     NA             
    ## [3,] "Lactobacillaceae" "Lactobacillus"
    ## [4,] "Lactobacillaceae" "Lactobacillus"
    ## [5,] "Lactobacillaceae" "Lactobacillus"
    ## [6,] "Lactobacillaceae" "Lactobacillus"

``` r
  taxTab
```

    ##                                                                                                                                                                                                                                                                                                                           Kingdom   
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Bacteria"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Bacteria"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Bacteria"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Bacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Bacteria"
    ##                                                                                                                                                                                                                                                                                                                           Phylum          
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         "Proteobacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   "Proteobacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Firmicutes"    
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Firmicutes"    
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Firmicutes"    
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Firmicutes"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Firmicutes"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Firmicutes"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Firmicutes"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Firmicutes"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Firmicutes"    
    ##                                                                                                                                                                                                                                                                                                                           Class                
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         "Alphaproteobacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   "Alphaproteobacteria"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Bacilli"            
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Bacilli"            
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Bacilli"            
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Bacilli"            
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Bacilli"            
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Bacilli"            
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Bacilli"            
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Bacilli"            
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Bacilli"            
    ##                                                                                                                                                                                                                                                                                                                           Order            
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         "Rickettsiales"  
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   "Rickettsiales"  
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillales"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Lactobacillales"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Lactobacillales"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Lactobacillales"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Lactobacillales"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Lactobacillales"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Lactobacillales"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Lactobacillales"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillales"
    ##                                                                                                                                                                                                                                                                                                                           Family            
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         "Mitochondria"    
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   "Mitochondria"    
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillaceae"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Lactobacillaceae"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Lactobacillaceae"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Lactobacillaceae"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Lactobacillaceae"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Lactobacillaceae"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Lactobacillaceae"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Lactobacillaceae"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillaceae"
    ##                                                                                                                                                                                                                                                                                                                           Genus          
    ## ATGAACTCCTACGGGAGGCAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACCCAAACTGGGATTAGAAACCCTTGTAGTCCAGAAGTG         NA             
    ## ATGAACTCCTACGGGAGGCAGCAGTGATAAACCTTTAGCAATAAACGAAAGTTTAACTGAGCTATACTAACTTTAGGGTTGGTTAATTTCATGCCAGCCACCACGGTCATACGATGAACCCAAGCTAATAGAGACTGGCGTAAAGAATGTTTTACATTATCCCTCAATAAAGCTAAATTTCACCTAAGTTGTAGAAAACCCTAGTTGATATAAAACAAACTACGAAAGTGGCTTTAATATTTCTGAATACACAATAGCGAAGATTCAAACTGGGATTAGAAACCCTTGTAGTCCT                   NA             
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillus"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTAGTAGTCCT          "Lactobacillus"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCGAGTAGTCCT          "Lactobacillus"
    ## TACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGAAACCCTTGTAGTCCT          "Lactobacillus"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGGGTAGTCCAGAAGTG "Lactobacillus"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCCTGTAGTCCAGAAGTG "Lactobacillus"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTTGTAGTCCAGAAGTG "Lactobacillus"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCGTGTAGTCCAGAAGTG "Lactobacillus"
    ## ATGAACTCCTACGGGAGGCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGTGCAGGCGGCTCGATAAGTCTGATGTGAAAGCCTTCGGCTCAACCGGAGAATTGCATCAGAAACTGTCGAGCTTGAGTACAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCGGTGGCGAAGGCGGCTCTCTGGTCTGTTACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGCATTAGAAACCCTAGTAGTCCAGAAGTG "Lactobacillus"

### Construct phylogenetic tree

``` r
  seqs <- getSequences(seqtabNoC) 
  names(seqs) <- seqs 
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

The phangorn R package is then used to construct a phylogenetic tree

``` r
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
     detach("package:phangorn", unload=TRUE)
```

``` r
  fitGTR
```

    ## 
    ##  loglikelihood: -1007.886 
    ## 
    ## unconstrained loglikelihood: -1052.16 
    ## Proportion of invariant sites: 0.1884369 
    ## Discrete gamma model
    ## Number of rate categories: 4 
    ## Shape parameter: 1.132012 
    ## 
    ## Rate matrix:
    ##              a            c         g         t
    ## a 0.000000e+00 1.108862e-05 0.8993070 0.8451276
    ## c 1.108862e-05 0.000000e+00 0.2137613 3.5146561
    ## g 8.993070e-01 2.137613e-01 0.0000000 1.0000000
    ## t 8.451276e-01 3.514656e+00 1.0000000 0.0000000
    ## 
    ## Base frequencies:  
    ## 0.3105401 0.2053405 0.2589221 0.2251972

``` r
  plot(fitGTR)
```

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Combine data into a phyloseq object

JE ne sais pas quoi mettre comme lien à la place de celui du workflow,
la suite de l’etude n’est alors pas bonne car elle ne prend pas les
données de l’article. J’ai quand même continué pour montrer mon
raisonnement.

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds") 
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

``` r
  samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE) 

  samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
  samdf <- samdf[!duplicated(samdf$SampleID),] # Remove duplicate entries for reverse reads
  rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy
  all(rownames(seqtabAll) %in% samdf$SampleID) # TRUE
```

    ## [1] FALSE

``` r
  rownames(samdf) <- samdf$SampleID
  keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment","host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass", "diet", "family_relationship", "genotype", "SampleID") 
  samdf <- samdf[rownames(seqtabAll), keep.cols]
```

``` r
#ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxTab), phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

### Loading the data

Le problème ici ets le même que celui mentionné precedemment. le
document ne correspond pas.

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

### Taxonomic filtering

``` r
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

Je ne n’ai pas les bonne données donc difficile de coparé mais si
j’avais eu les données j’aurais comparé combien d’ASV ont été trouvé
comparé aux OTU trouvés dans l’article. Dans l’article, 138 OTU ont été
identifiées dans les deux groupes d’étude.

``` r
prevdf = apply(X = otu_table(ps),MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps),tax_table(ps))
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

## Diversité alpha

Dans cette partie je dois calculer les especes observé et les 2 indices
de l’alpha diversité: L’indice de Shannon qui permet de regarder la
diverité des population du microbiote( il est inversement proportionel à
la divesité) et l’indice de Chao 1 qui est un extimateur de richesse.

Il existe different packages qui permettent de mesurere les indices
d’alpha diversité tel que le package vegan, mais je vais utiliser
phyloseq et ggplot2 pour la mise en page des graphiques.

avant de proceder à faire des plots d’alpha diversité il faut creer un
nouvel objet qui contient les données de ps mais de manière à ce que
tous les echantillions ne soit pas sur graph mais qu’il y ait les
quantiles pour l’ensemble des echantillions. Je n’ai deja as les bonnes
données vu que j’ai celle du wrokflow. Aussi je ne saurais pas comment à
partir de l’objet ps prendre les 30 premier qui correspondent aux femmes
qui n’ont pas d’infection par le HPV et les 30 autres en “early stage of
infection”. ETanta données qu’il ya. des erreurs dans mon code je ne
peux pas faire la suite de ce que je voulais faire. J’aurais aimé
pouvoir refaire cette analyse en faisant cette fois avec la methode
d’OTU en plus, et ainsi à l’aide de wilcox.test() comparer les 3 indices
des 2 groupes (HPV- et HPV+) et entre les 2 methodes en regardant si la
p-value est superieur ou inferieur au threshold à 5%.

``` r
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 383 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 383 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 383 tips and 381 internal nodes ]

``` r
estimate_richness(ps, split = TRUE, measures = c("Observed", "Chao1", "Shannon"))
```

    ##        Observed     Chao1  se.chao1   Shannon
    ## F3D0        222 256.44000 13.989381 4.0731879
    ## F3D1        202 232.87500 13.010810 4.2171130
    ## F3D11       209 268.23077 21.073720 2.9264568
    ## F3D125      258 286.92308 11.000894 3.9238666
    ## F3D13       185 229.20000 16.135918 2.9741181
    ## F3D141      186 255.39130 24.769718 3.6321593
    ## F3D142      159 229.03704 23.803228 3.4741998
    ## F3D143      155 208.04000 19.536828 3.6288149
    ## F3D144      168 233.03226 21.576653 3.3666345
    ## F3D145      193 246.20000 18.615647 3.3033210
    ## F3D146      213 279.12000 23.254405 3.9427854
    ## F3D147      238 266.75000 11.160002 3.4542081
    ## F3D148      226 279.00000 19.311814 3.5839469
    ## F3D149      227 264.62500 15.131786 3.7805745
    ## F3D15       182 218.75000 13.790513 3.1386157
    ## F3D150      206 261.12000 20.136450 3.9085238
    ## F3D17       230 273.36364 15.509094 3.5413999
    ## F3D19       159 227.87500 24.313458 3.5089232
    ## F3D2        243 269.25000 11.511631 3.6194402
    ## F3D21       260 299.60000 15.568389 3.9868846
    ## F3D25       174 247.00000 22.648357 3.2777465
    ## F3D3        161 188.36364 10.985242 3.2047233
    ## F3D364      184 247.17143 20.390461 3.2256392
    ## F3D5        194 243.00000 18.362526 4.0544366
    ## F3D6        223 273.09091 17.317321 3.6255099
    ## F3D65       141 209.87500 24.313329 3.2558176
    ## F3D7        163 209.16129 16.537500 3.1943057
    ## F3D8        217 252.68966 13.831752 3.9864473
    ## F3D9        225 246.96875  9.421142 3.9946995
    ## F4D0        210 253.96552 16.215765 4.0707294
    ## F4D1        202 227.10714 10.722193 4.1815938
    ## F4D11       163 260.65000 34.141788 2.9213434
    ## F4D125      264 284.18182  8.798310 3.7578242
    ## F4D13       184 231.00000 17.977473 3.3002048
    ## F4D141      249 269.58333  8.753587 3.8496207
    ## F4D142      256 275.19355  8.604540 3.9988990
    ## F4D143      286 296.00000  5.877400 3.9204484
    ## F4D144      270 313.38462 16.528133 3.6672217
    ## F4D145      280 291.66667  6.785391 4.0082577
    ## F4D146      118 190.05882 27.987115 2.9395011
    ## F4D147      257 271.50000  7.090803 4.0017242
    ## F4D148      268 291.15625  9.789587 3.9396166
    ## F4D149      229 265.02941 13.381104 3.3677896
    ## F4D15       185 238.04000 19.536983 3.3784129
    ## F4D150      283 310.55556 12.859898 3.9240880
    ## F4D17       165 202.12121 13.787509 3.2981457
    ## F4D19       185 218.60976 12.140055 3.5160816
    ## F4D2        151 196.00000 17.583875 3.2549000
    ## F4D21       133 184.13043 19.433968 2.8822837
    ## F4D25       193 261.44000 23.902062 3.1220057
    ## F4D3        173 268.05556 34.527560 3.1671616
    ## F4D302      296 314.60000  8.832856 4.1217008
    ## F4D4        222 248.71429 12.065851 3.6122604
    ## F4D5        188 227.80769 15.468067 3.1830155
    ## F4D6        186 229.00000 17.181671 3.2009632
    ## F4D65       157 218.60000 21.982783 3.2161943
    ## F4D7        177 246.78947 26.354923 3.1051385
    ## F4D8        206 243.84000 15.033780 3.8517740
    ## F4D9        201 227.10526 12.166313 4.1025018
    ## F5D0        238 268.87500 13.010908 3.8466073
    ## F5D1        222 257.03704 13.892878 3.5463776
    ## F5D11       229 252.78571 10.302176 3.1015775
    ## F5D125      222 267.37037 16.935657 4.1912881
    ## F5D13       251 271.51724  9.168474 3.1294512
    ## F5D141      190 281.50000 32.361054 3.8542358
    ## F5D142      164 257.88235 34.803179 3.4911106
    ## F5D143      266 291.14286 11.525537 3.4746850
    ## F5D144      220 298.71429 28.186109 3.9449704
    ## F5D145       32  58.25000 18.743970 2.7547682
    ## F5D146      215 272.18750 19.344475 3.9574594
    ## F5D147      124 154.37037 12.469707 3.2513771
    ## F5D148      177 232.31250 18.850443 3.4468610
    ## F5D149      178 231.00000 19.109273 3.3039562
    ## F5D15       159 228.78947 26.354769 3.2593122
    ## F5D150      192 263.86957 25.477597 3.8490144
    ## F5D165       22  44.00000 17.432989 0.2994423
    ## F5D17       210 240.00000 12.473616 3.7667466
    ## F5D19       149 198.00000 18.362302 3.3370754
    ## F5D2        224 251.39130 12.014251 3.7256759
    ## F5D21       104 174.71429 29.270741 2.8198628
    ## F5D25       188 230.50000 15.658648 3.8344834
    ## F5D3        141 194.26087 20.069059 2.8004808
    ## F5D364      149 166.65217  8.667588 3.2728420
    ## F5D4        177 216.13636 15.969012 3.2207894
    ## F5D45       167 193.27778 10.447960 3.3576673
    ## F5D5        221 240.89474  9.913286 4.0573178
    ## F5D6        196 216.66667  9.635127 3.0106320
    ## F5D65       132 179.04000 17.787094 3.0406259
    ## F5D7        230 265.65217 14.683947 3.4482826
    ## F5D8        205 241.03333 13.812019 3.3771702
    ## F5D9        242 282.28571 15.306757 3.4120079
    ## F6D0        192 223.13793 12.480249 3.6921785
    ## F6D1        254 278.79167 11.029376 4.1794408
    ## F6D11       164 222.80000 22.607547 2.8063342
    ## F6D125      223 294.29167 24.994436 3.4064385
    ## F6D13       188 229.40000 16.110789 3.2641215
    ## F6D141      236 270.73077 13.936178 3.7134966
    ## F6D142      183 232.00000 18.362479 3.1037035
    ## F6D143      179 218.51613 14.687667 3.0427898
    ## F6D144      237 254.53125  8.005827 3.7926820
    ## F6D145      137 186.13636 19.073663 3.1649615
    ## F6D146      211 258.11538 17.618527 3.3306169
    ## F6D147      249 313.47368 24.719219 4.0677277
    ## F6D148      198 234.38710 13.798500 4.0164536
    ## F6D149      204 259.00000 19.472063 3.3657108
    ## F6D15       144 187.93750 19.123568 3.1184102
    ## F6D150      162 205.55556 16.410123 3.5580041
    ## F6D165       39  69.00000 17.883682 0.2701448
    ## F6D17       169 220.00000 18.739536 3.0913877
    ## F6D19       223 266.12500 16.813400 3.7602833
    ## F6D2        220 263.38462 16.528010 3.7201310
    ## F6D21       201 244.06250 15.551175 3.1164215
    ## F6D25       256 284.50000 12.001720 4.0262079
    ## F6D3        257 267.68750  6.540310 4.1737206
    ## F6D364       85 103.90000  9.422257 2.8265433
    ## F6D4        197 230.17647 12.580543 3.2378632
    ## F6D45       204 251.22222 17.468180 3.4255232
    ## F6D5        221 258.05000 15.708016 4.1941735
    ## F6D6        167 231.56522 23.381153 3.1754119
    ## F6D65       255 290.41667 13.024546 3.6458281
    ## F6D7        171 202.13793 12.480188 3.0731338
    ## F6D8        198 255.11538 20.478733 3.4902785
    ## F6D9        231 266.15000 15.079905 3.5060535
    ## F7D0        209 231.89474 11.015143 3.9252533
    ## F7D1        221 243.88462 10.191560 3.8461478
    ## F7D11       203 230.03846 11.537388 3.5484879
    ## F7D125      232 252.21739  9.577071 3.7690177
    ## F7D13       206 241.45455 14.795771 3.6052391
    ## F7D141      175 204.28571 12.024356 3.5446711
    ## F7D142      146 199.04000 19.536772 3.2811761
    ## F7D143      169 226.03704 20.239506 3.0297847
    ## F7D144      183 230.51724 17.215165 3.1718134
    ## F7D145      200 266.12000 23.254350 3.7533233
    ## F7D146      210 233.78571 10.302140 3.7768157
    ## F7D147      176 239.75000 24.121915 3.6173382
    ## F7D148      220 257.84000 15.033822 3.4473971
    ## F7D149      147 209.21739 22.700422 3.0788379
    ## F7D15       169 210.13043 16.396304 2.8323945
    ## F7D150      159 186.03846 11.537248 3.2251028
    ## F7D2        225 256.95455 13.661821 3.8458921
    ## F7D25       133 198.33333 25.395095 2.7242756
    ## F7D3        239 261.14286 10.476563 4.2159215
    ## F7D4        173 195.03704  9.822997 3.5708744
    ## F7D45       181 251.83333 27.118341 3.4323908
    ## F7D5        143 167.23077 10.632135 3.2417468
    ## F7D6        209 231.28571  9.327969 3.6577014
    ## F7D65       221 266.72414 16.712268 3.2505425
    ## F7D7        178 218.55172 15.242365 3.3032809
    ## F7D9        185 212.33333 11.230202 3.3535029
    ## F8D0        174 225.25000 21.627694 3.5700491
    ## F8D1        163 198.45455 14.795578 3.1051087
    ## F8D125      210 259.28571 19.373871 3.8650241
    ## F8D141      140 174.44000 13.989059 3.3547359
    ## F8D142      173 236.37037 21.989028 3.3172239
    ## F8D143      190 227.27586 14.295348 3.6058112
    ## F8D144      189 253.47368 24.718881 3.9331624
    ## F8D145      155 206.33333 18.107142 3.2338450
    ## F8D146      188 228.88571 14.609535 3.6426245
    ## F8D147      189 229.28571 15.306616 3.6076570
    ## F8D148      202 262.27273 22.422938 3.8176744
    ## F8D149      126 163.00000 15.913309 3.5094185
    ## F8D150      204 230.40000 12.108933 3.8817638
    ## F8D2        153 209.40000 21.867063 3.0770867
    ## F8D25       188 270.50000 30.719523 3.5200987
    ## F8D3        157 192.25000 13.364525 2.9800151
    ## F8D4        150 217.56250 27.052756 3.0566958
    ## F8D45       161 204.04348 16.985510 3.5950396
    ## F8D5        187 211.23077 10.632289 3.8050391
    ## F8D6        123 172.50000 19.712683 2.8903128
    ## F8D65       208 229.96875  9.421116 3.8206348
    ## F8D7        171 293.50000 51.936945 3.5358283
    ## F8D8        104 147.05000 17.659229 2.5418595
    ## F8D9        182 221.51613 14.687676 3.0359410
    ## M1D0        221 238.71429  8.307524 3.5983778
    ## M1D1        160 188.33333 12.616635 3.4150778
    ## M1D11       184 249.80769 22.904119 2.8022951
    ## M1D125      200 246.94118 16.347064 3.4075460
    ## M1D13       210 238.33333 12.616845 3.5963395
    ## M1D141      134 213.68750 30.977716 3.7332683
    ## M1D142      210 237.02500 10.408490 3.5023118
    ## M1D143      134 204.71429 25.842355 2.8749473
    ## M1D144      232 262.51613 12.091892 3.7056064
    ## M1D145      238 269.53846 12.953437 3.6975157
    ## M1D146      150 218.87500 24.313396 3.4445996
    ## M1D147      171 212.62162 14.600003 3.1106106
    ## M1D148      218 251.78125 12.944483 3.6738034
    ## M1D149       18  27.00000  8.025350 2.6984822
    ## M1D15       170 223.04000 19.536910 3.3034635
    ## M1D17       114 175.50000 26.100177 2.9930800
    ## M1D19       168 246.15789 28.898965 3.3475898
    ## M1D2        123 160.14286 15.530049 3.1960309
    ## M1D21       132 177.00000 17.802318 3.1920362
    ## M1D25       200 243.96552 16.215736 3.4714127
    ## M1D3        213 253.61538 19.076722 3.5839957
    ## M1D364      175 221.24324 15.809715 3.3260573
    ## M1D4        193 218.58824 12.334873 3.0841705
    ## M1D5        106 163.40000 24.160438 3.0933310
    ## M1D6        191 207.91667  8.320458 3.3236833
    ## M1D65       171 224.26087 20.069258 3.7185776
    ## M1D7        113 154.16667 17.562207 3.3374783
    ## M1D8        196 229.05556 14.796248 3.2022838
    ## M1D9         17  44.50000 22.663451 2.6318526
    ## M2D0        200 239.13636 15.969116 4.1379927
    ## M2D1        180 211.95455 13.661659 3.6999095
    ## M2D11       124 196.76923 30.678297 3.1998478
    ## M2D125       61 110.58333 22.943413 3.2451120
    ## M2D13       132 202.50000 28.010716 2.7433081
    ## M2D141      171 218.00000 17.977410 3.4814903
    ## M2D142      206 251.53571 16.817114 3.4482954
    ## M2D143      146 181.65217 14.683594 3.2085748
    ## M2D144      205 235.44118 11.801022 3.3606277
    ## M2D145      207 298.63636 31.449372 3.4808099
    ## M2D146      187 219.34375 12.529933 3.4757944
    ## M2D147      199 237.07692 14.949722 3.4933818
    ## M2D148      175 212.27586 14.295300 3.3311760
    ## M2D149      158 194.03333 13.811876 3.1263696
    ## M2D15       151 192.13043 16.396201 2.8581031
    ## M2D150      191 226.33333 12.749734 3.1867648
    ## M2D17       132 179.04545 18.433086 2.7038866
    ## M2D19        18  34.50000 12.877005 2.7502140
    ## M2D2        117 168.75000 20.419450 2.9114302
    ## M2D21       177 266.37500 29.997157 3.0266665
    ## M2D25       201 235.02778 12.641905 3.2675560
    ## M2D3        196 209.04348  6.962774 3.6366126
    ## M2D364      251 281.00000 12.473706 3.8646979
    ## M2D4        226 248.00000 10.091102 3.5743542
    ## M2D5        207 232.20000 11.052741 3.0634564
    ## M2D6        235 260.50000 11.512374 3.5879035
    ## M2D65       256 312.43750 23.374931 3.8365409
    ## M2D7        202 232.51613 12.091834 3.3500365
    ## M2D8        208 234.89655 11.188394 3.5644380
    ## M2D9        217 241.24138 10.360379 3.3946096
    ## M3D0        205 231.64000 11.521597 3.4512232
    ## M3D1        142 204.66667 24.553305 3.4131892
    ## M3D11       213 259.66667 16.402351 3.5308988
    ## M3D125      205 247.00000 15.803252 3.2900242
    ## M3D13       227 240.03448  6.623122 3.3567222
    ## M3D141      178 205.02857 10.736007 3.7714875
    ## M3D142      226 247.02564  8.733855 3.3524701
    ## M3D143      249 275.09091 10.607055 3.4065596
    ## M3D144      215 260.93333 16.619513 3.7024050
    ## M3D145      235 277.24138 15.725898 3.5810235
    ## M3D146      185 244.42857 19.446659 3.2948363
    ## M3D147      171 202.53846 12.953247 3.9628005
    ## M3D148      126 213.35294 32.785108 3.8360241
    ## M3D149        2   2.00000  0.000000 0.6931472
    ## M3D15       165 189.23077 10.632221 3.5776812
    ## M3D150      197 254.00000 19.828107 3.4812469
    ## M3D17       215 283.07692 23.529621 3.4086795
    ## M3D175       89  94.00000  4.132693 2.9707890
    ## M3D19       193 256.57692 22.286387 3.3047918
    ## M3D2        218 249.07143 15.108387 4.0303024
    ## M3D21       195 233.60714 14.816955 2.9963481
    ## M3D25       175 210.68966 13.831636 3.3668902
    ## M3D3         95 109.13043  7.374338 3.6110432
    ## M3D364      162 240.00000 35.344870 3.6906515
    ## M3D4        193 252.40000 21.359381 3.7755953
    ## M3D45       189 234.60000 15.869187 3.5751432
    ## M3D5         53  84.23077 15.482908 3.5768338
    ## M3D6        201 260.40000 21.359415 3.4774218
    ## M3D65       260 290.37037 12.470139 3.9996313
    ## M3D7        168 194.55882 10.671086 3.1224019
    ## M3D8        109 176.57143 28.195223 3.5968707
    ## M3D9        177 234.11538 20.478642 3.4749801
    ## M4D0        164 184.81250  9.058450 3.4724661
    ## M4D1        145 176.88889 12.936752 3.2475506
    ## M4D11       170 211.00000 16.769565 3.7085549
    ## M4D125      250 284.16667 14.054028 4.1575711
    ## M4D13       179 232.20000 18.615599 3.2016005
    ## M4D141      195 276.47619 28.987705 3.5863120
    ## M4D142      250 290.55172 15.242552 3.7278672
    ## M4D143      259 275.46667  7.134638 3.8644118
    ## M4D144      244 271.51220 10.485396 4.0048369
    ## M4D145      186 228.00000 15.803192 3.9438196
    ## M4D146      221 286.20690 22.035548 3.8049485
    ## M4D147      214 254.52941 14.621005 3.9519968
    ## M4D148      151 194.67647 15.473486 3.4722616
    ## M4D149      177 248.50000 23.469077 3.3970922
    ## M4D15       181 236.12000 20.136343 3.2890188
    ## M4D150      173 222.87500 17.402510 3.7017665
    ## M4D17       173 218.37037 16.935489 3.3667366
    ## M4D175      102 119.10000 10.410811 2.5876471
    ## M4D19       167 198.53846 12.953231 3.3765158
    ## M4D2        171 205.16667 14.053778 3.8162313
    ## M4D21       194 224.93750 12.121196 3.7905852
    ## M4D25       182 241.00000 20.370155 3.8069273
    ## M4D3        222 250.88889 12.010469 4.0088640
    ## M4D364      169 199.56522 13.054148 3.8404872
    ## M4D4        135 182.11538 17.618174 3.0617213
    ## M4D45       197 275.40000 31.192139 4.2100174
    ## M4D5        207 244.27586 14.295395 4.0042297
    ## M4D6        132 171.26087 15.816011 3.1524843
    ## M4D65       190 239.04348 18.808429 3.8036067
    ## M4D7        190 243.12500 19.788352 3.5514237
    ## M4D8        156 187.88889 12.936801 3.5808184
    ## M4D9        179 202.90323 10.093754 3.7533304
    ## M5D0        179 207.50000 12.001534 3.7746287
    ## M5D1        171 216.00000 17.583992 3.2723046
    ## M5D11       234 294.88235 24.403557 3.9294051
    ## M5D125      235 269.16667 14.053991 3.9267102
    ## M5D13       212 232.32258  8.967678 3.5853162
    ## M5D141      218 258.83333 15.187686 3.4679273
    ## M5D142      199 231.66667 12.264165 3.2835319
    ## M5D143      220 262.24138 15.725863 3.9339038
    ## M5D144      184 221.27586 14.295330 3.6850371
    ## M5D145      232 262.27273 13.109628 3.9792529
    ## M5D146      233 251.03030  8.113727 3.7415212
    ## M5D147      233 254.79412  9.238736 3.6183955
    ## M5D148      264 290.03704 11.113433 4.0006868
    ## M5D149      205 282.53846 26.106939 3.5235416
    ## M5D15       225 262.00000 15.913942 4.1870955
    ## M5D150      239 286.04000 17.787575 3.8872187
    ## M5D17       203 240.84000 15.033770 3.4959084
    ## M5D175       65  80.16667 10.854479 2.9861585
    ## M5D19       223 256.38710 12.933233 3.1555672
    ## M5D2        236 251.00000  7.949566 4.1281358
    ## M5D21       269 290.56522 10.045871 4.0968201
    ## M5D25       155 187.62069 12.924006 3.0732802
    ## M5D3        225 279.05000 21.138286 3.8388832
    ## M5D364      178 187.24000  5.365277 3.9138368
    ## M5D4        184 250.12000 23.254276 3.4383836
    ## M5D45       214 236.20000  9.635501 3.8871725
    ## M5D5        214 249.45455 14.795800 3.4784802
    ## M5D6        185 248.75000 24.121977 3.6647265
    ## M5D65       195 226.36364 12.153409 3.0197664
    ## M5D7        211 274.37037 21.989186 3.5855551
    ## M5D8        221 244.80000 10.592174 4.0256820
    ## M5D9        204 242.89655 14.765632 3.6736996
    ## M6D0        223 268.88235 19.450665 3.8868199
    ## M6D1        214 251.60000 14.264328 3.8364131
    ## M6D11       205 250.12000 17.220474 4.0316933
    ## M6D125      242 271.29167 12.502426 3.7115903
    ## M6D13       142 166.24138 10.360178 3.0320918
    ## M6D141      209 239.15385 11.342015 3.4733326
    ## M6D142      193 242.11111 18.007967 3.5104394
    ## M6D143      185 253.90000 25.680849 4.0030362
    ## M6D144      198 276.71429 28.185985 3.7333670
    ## M6D145      232 264.80000 13.479447 4.0576663
    ## M6D146      229 249.77778  9.407874 4.2376865
    ## M6D147      224 258.50000 13.366042 3.9710225
    ## M6D148      201 231.44118 11.801014 3.9443322
    ## M6D149      214 261.14286 18.707175 3.9998170
    ## M6D15       168 202.50000 13.365895 3.0625102
    ## M6D150      229 280.10714 18.389566 3.5469472
    ## M6D17       167 222.10000 19.130108 3.2066879
    ## M6D175       86 105.00000 11.612049 3.2234678
    ## M6D19       141 215.25000 27.283871 2.8679215
    ## M6D2        215 229.13043  7.374791 3.6602844
    ## M6D21       173 200.33333 11.230170 3.0194063
    ## M6D25       181 230.03846 18.175178 3.6204846
    ## M6D3        197 257.88235 24.403315 3.4806403
    ## M6D364      184 197.12500  9.011669 4.0816790
    ## M6D4        196 239.00000 17.181717 2.8972594
    ## M6D45       214 248.73077 13.936123 3.8049882
    ## M6D5        238 252.66667  6.882462 3.5113829
    ## M6D6        173 230.03333 19.650882 3.0823013
    ## M6D65       251 283.50000 13.528136 3.6394325
    ## M6D7        217 305.66667 32.597980 2.9135561
    ## M6D8        210 265.12000 20.136465 3.7814708
    ## M6D9        160 214.47368 21.590587 3.0819607

``` r
plot_richness(ps, measures=c("Observed", "Chao1", "Shannon"), title = " Alpha diversity indices for the microbiotes from both non-infected women and women in early stageof infection my the Human Papillomavirus.")
```

![](CC2_HPV_traitement_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->
