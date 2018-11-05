Analysis Report 1: Amplicons and the Letter Q
================
Derek Newberger
November 11, 2018

Introduction
============

Add about 1.5-2 pages here. Must cite at least 5 peer reviewed articles.

Methods
=======

Sample origin and sequencing
----------------------------

Add about half a page here. In this section instead of first person (I/we), use Fierer et al., since you'll just be describing what they did, based on the methods in their paper.

Computational
-------------

These are the methods you used. Should probably be at least a half of a page. At a very minimum should include citations for DADA2 (Callahan *et al.*, 2016) and phyloseq (McMurdie and Holmes, 2013). Note that these don't count towards the five references you need to cite in the introduction. We will use the Ribosomal Database Project reference taxonomy to assign taxonomy to unkown sequences (Cole *et al.*, 2008).

Results
=======

In addition to a minimum of 3-4 figures/tables (and associated captions), you should include sufficient text in this section to describe what your findings were. Remember that in the results section you just describe what you found, but you don't interpret it - that happens in the discussion.

``` r
# Be sure to install these packages before running this script
# They can be installed either with the intall.packages() function
# or with the 'Packages' pane in RStudio

# load general-use packages
library("dplyr")
library("tidyr")
library("knitr")
library("ggplot2")

# this package allows for the easy inclusion of literature citations in our Rmd
# more info here: https://github.com/crsh/citr
# and here:
# http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html
library("citr")

# These are the primary packages well use to clean and analyze the data
# this package needs to be installed from bioconductor -- it's not on CRAN
# see info here: https://benjjneb.github.io/dada2/dada-installation.html
library("dada2")

# This to export a fasta of our final denoised sequence variants
library("seqinr")

# To install this you have to install from GitHub
# See more info here: https://github.com/leffj/mctoolsr
# run this -- install.packages("devtools")
# and then this -- devtools::install_github("leffj/mctoolsr")
library("mctoolsr")

# And this to visualize our results
# it also needs to be installed from bioconductor
library("phyloseq")
```

``` r
# NOTE: Much of the following follows the DADA2 tutorials available here:
# https://benjjneb.github.io/dada2/tutorial.html
# Accessed October 19, 2017

# set the base path for our input data files
path <- "data/raw_data"

# Sort ensures samples are in order
filenames_forward_reads <- sort(list.files(path, pattern = ".fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME.fastq
sample_names <- sapply(strsplit(filenames_forward_reads, "\\."), `[`, 1)

# Specify the full path to each of the filenames_forward_reads
filenames_forward_reads <- file.path(path, filenames_forward_reads)
```

We can see from the quality profiles that most reads tend to get pretty bad in quality after around 200 bases.

``` r
# Place filtered files in filtered/ subdirectory
# note this will fail if the directory doesn't exist
filter_path <- file.path("output", "filtered")
filtered_reads_path <- file.path(filter_path,
                                 paste0(sample_names,
                                        "_filt.fastq.gz"))

# See ?filterAndTrim for details on the parameters
# See here for adjustments for 454 data:
# https://benjjneb.github.io/dada2/
#     faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
filtered_output <- filterAndTrim(fwd = filenames_forward_reads,
                                 filt = filtered_reads_path,
                                 maxLen = 300,
                                 maxN = 0, # discard any seqs with Ns
                                 maxEE = 3, # allow w/ up to 3 expected errors
                                 truncQ = 2, # cut off if quality gets this low
                                 rm.phix = TRUE,
                                 compress = TRUE,
                                 multithread = FALSE)
```

``` r
# this build error models from each of the samples
errors_forward_reads <- learnErrors(filtered_reads_path,
                                   multithread = FALSE)
```

    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## Not all sequences were the same length.
    ## 28802208 total bases in 125531 reads from 115 samples will be used for learning the error rates.
    ## Initializing error rates to maximum possible estimate.
    ## selfConsist step 1 ...................................................................................................................
    ##    selfConsist step 2
    ##    selfConsist step 3
    ##    selfConsist step 4
    ##    selfConsist step 5
    ## Convergence after  5  rounds.

``` r
# get rid of any duplicated sequences
dereplicated_forward_reads <- derepFastq(filtered_reads_path,
                                         verbose = TRUE)
```

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938262_filt.fastq.gz

    ## Encountered 139 unique sequences from 318 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938263_filt.fastq.gz

    ## Encountered 416 unique sequences from 895 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938264_filt.fastq.gz

    ## Encountered 317 unique sequences from 663 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938265_filt.fastq.gz

    ## Encountered 117 unique sequences from 216 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938266_filt.fastq.gz

    ## Encountered 135 unique sequences from 325 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938267_filt.fastq.gz

    ## Encountered 166 unique sequences from 563 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938268_filt.fastq.gz

    ## Encountered 212 unique sequences from 571 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938269_filt.fastq.gz

    ## Encountered 607 unique sequences from 1267 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938270_filt.fastq.gz

    ## Encountered 704 unique sequences from 1396 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938271_filt.fastq.gz

    ## Encountered 697 unique sequences from 1204 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938272_filt.fastq.gz

    ## Encountered 614 unique sequences from 1034 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938273_filt.fastq.gz

    ## Encountered 662 unique sequences from 1079 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938274_filt.fastq.gz

    ## Encountered 818 unique sequences from 1543 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938275_filt.fastq.gz

    ## Encountered 744 unique sequences from 1239 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938276_filt.fastq.gz

    ## Encountered 696 unique sequences from 1352 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938277_filt.fastq.gz

    ## Encountered 592 unique sequences from 1162 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938278_filt.fastq.gz

    ## Encountered 749 unique sequences from 1533 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938279_filt.fastq.gz

    ## Encountered 528 unique sequences from 1102 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938280_filt.fastq.gz

    ## Encountered 530 unique sequences from 956 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938281_filt.fastq.gz

    ## Encountered 580 unique sequences from 1124 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938282_filt.fastq.gz

    ## Encountered 523 unique sequences from 1099 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938283_filt.fastq.gz

    ## Encountered 638 unique sequences from 1355 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938284_filt.fastq.gz

    ## Encountered 600 unique sequences from 1263 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938285_filt.fastq.gz

    ## Encountered 536 unique sequences from 1261 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938286_filt.fastq.gz

    ## Encountered 659 unique sequences from 1308 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938287_filt.fastq.gz

    ## Encountered 492 unique sequences from 1106 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938288_filt.fastq.gz

    ## Encountered 251 unique sequences from 474 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938289_filt.fastq.gz

    ## Encountered 551 unique sequences from 1240 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938290_filt.fastq.gz

    ## Encountered 440 unique sequences from 970 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938291_filt.fastq.gz

    ## Encountered 579 unique sequences from 1215 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938292_filt.fastq.gz

    ## Encountered 401 unique sequences from 691 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938293_filt.fastq.gz

    ## Encountered 486 unique sequences from 820 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938294_filt.fastq.gz

    ## Encountered 419 unique sequences from 1160 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938295_filt.fastq.gz

    ## Encountered 505 unique sequences from 1059 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938296_filt.fastq.gz

    ## Encountered 331 unique sequences from 620 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938297_filt.fastq.gz

    ## Encountered 668 unique sequences from 1215 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938298_filt.fastq.gz

    ## Encountered 640 unique sequences from 1450 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938299_filt.fastq.gz

    ## Encountered 755 unique sequences from 1219 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938300_filt.fastq.gz

    ## Encountered 531 unique sequences from 1232 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938301_filt.fastq.gz

    ## Encountered 515 unique sequences from 943 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938302_filt.fastq.gz

    ## Encountered 590 unique sequences from 1392 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938303_filt.fastq.gz

    ## Encountered 504 unique sequences from 1065 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938304_filt.fastq.gz

    ## Encountered 633 unique sequences from 1230 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938305_filt.fastq.gz

    ## Encountered 554 unique sequences from 1049 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938306_filt.fastq.gz

    ## Encountered 296 unique sequences from 502 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938307_filt.fastq.gz

    ## Encountered 667 unique sequences from 1384 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938308_filt.fastq.gz

    ## Encountered 625 unique sequences from 1112 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938309_filt.fastq.gz

    ## Encountered 227 unique sequences from 674 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938310_filt.fastq.gz

    ## Encountered 291 unique sequences from 1125 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938311_filt.fastq.gz

    ## Encountered 238 unique sequences from 1023 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938312_filt.fastq.gz

    ## Encountered 309 unique sequences from 1189 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938313_filt.fastq.gz

    ## Encountered 385 unique sequences from 1273 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938314_filt.fastq.gz

    ## Encountered 237 unique sequences from 809 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938315_filt.fastq.gz

    ## Encountered 342 unique sequences from 1098 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938316_filt.fastq.gz

    ## Encountered 453 unique sequences from 1234 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938317_filt.fastq.gz

    ## Encountered 185 unique sequences from 529 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938318_filt.fastq.gz

    ## Encountered 88 unique sequences from 241 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938319_filt.fastq.gz

    ## Encountered 440 unique sequences from 1328 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938320_filt.fastq.gz

    ## Encountered 292 unique sequences from 723 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938321_filt.fastq.gz

    ## Encountered 543 unique sequences from 1947 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938322_filt.fastq.gz

    ## Encountered 527 unique sequences from 1416 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938323_filt.fastq.gz

    ## Encountered 363 unique sequences from 1127 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938324_filt.fastq.gz

    ## Encountered 225 unique sequences from 752 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938325_filt.fastq.gz

    ## Encountered 530 unique sequences from 1320 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938326_filt.fastq.gz

    ## Encountered 618 unique sequences from 1685 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938327_filt.fastq.gz

    ## Encountered 251 unique sequences from 753 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938328_filt.fastq.gz

    ## Encountered 362 unique sequences from 1022 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938329_filt.fastq.gz

    ## Encountered 277 unique sequences from 1245 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938330_filt.fastq.gz

    ## Encountered 598 unique sequences from 1531 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938331_filt.fastq.gz

    ## Encountered 426 unique sequences from 1160 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938332_filt.fastq.gz

    ## Encountered 395 unique sequences from 1143 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938333_filt.fastq.gz

    ## Encountered 369 unique sequences from 1211 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938334_filt.fastq.gz

    ## Encountered 249 unique sequences from 831 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938335_filt.fastq.gz

    ## Encountered 470 unique sequences from 1489 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938336_filt.fastq.gz

    ## Encountered 374 unique sequences from 1409 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938337_filt.fastq.gz

    ## Encountered 280 unique sequences from 1075 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938338_filt.fastq.gz

    ## Encountered 332 unique sequences from 1332 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938339_filt.fastq.gz

    ## Encountered 360 unique sequences from 1127 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938340_filt.fastq.gz

    ## Encountered 386 unique sequences from 1258 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938341_filt.fastq.gz

    ## Encountered 165 unique sequences from 499 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938342_filt.fastq.gz

    ## Encountered 256 unique sequences from 1184 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938343_filt.fastq.gz

    ## Encountered 310 unique sequences from 1017 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938344_filt.fastq.gz

    ## Encountered 273 unique sequences from 1109 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938345_filt.fastq.gz

    ## Encountered 317 unique sequences from 1285 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938346_filt.fastq.gz

    ## Encountered 294 unique sequences from 1067 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938347_filt.fastq.gz

    ## Encountered 348 unique sequences from 1049 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938348_filt.fastq.gz

    ## Encountered 289 unique sequences from 1056 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938349_filt.fastq.gz

    ## Encountered 358 unique sequences from 1493 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938350_filt.fastq.gz

    ## Encountered 336 unique sequences from 1348 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938351_filt.fastq.gz

    ## Encountered 238 unique sequences from 742 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938352_filt.fastq.gz

    ## Encountered 482 unique sequences from 2061 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938353_filt.fastq.gz

    ## Encountered 142 unique sequences from 452 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938354_filt.fastq.gz

    ## Encountered 397 unique sequences from 1499 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938355_filt.fastq.gz

    ## Encountered 524 unique sequences from 2281 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938356_filt.fastq.gz

    ## Encountered 389 unique sequences from 1338 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938357_filt.fastq.gz

    ## Encountered 277 unique sequences from 1090 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938358_filt.fastq.gz

    ## Encountered 239 unique sequences from 902 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938359_filt.fastq.gz

    ## Encountered 489 unique sequences from 1584 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938360_filt.fastq.gz

    ## Encountered 328 unique sequences from 1052 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938361_filt.fastq.gz

    ## Encountered 322 unique sequences from 1277 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938362_filt.fastq.gz

    ## Encountered 130 unique sequences from 609 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938363_filt.fastq.gz

    ## Encountered 509 unique sequences from 2013 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938364_filt.fastq.gz

    ## Encountered 353 unique sequences from 1252 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938365_filt.fastq.gz

    ## Encountered 187 unique sequences from 819 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938366_filt.fastq.gz

    ## Encountered 92 unique sequences from 277 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938367_filt.fastq.gz

    ## Encountered 327 unique sequences from 1077 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938368_filt.fastq.gz

    ## Encountered 577 unique sequences from 1745 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938369_filt.fastq.gz

    ## Encountered 378 unique sequences from 1528 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938370_filt.fastq.gz

    ## Encountered 260 unique sequences from 1191 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938371_filt.fastq.gz

    ## Encountered 277 unique sequences from 996 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938372_filt.fastq.gz

    ## Encountered 243 unique sequences from 841 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938373_filt.fastq.gz

    ## Encountered 377 unique sequences from 982 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938374_filt.fastq.gz

    ## Encountered 272 unique sequences from 462 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938375_filt.fastq.gz

    ## Encountered 267 unique sequences from 718 total sequences read.

    ## Not all sequences were the same length.

    ## Dereplicating sequence entries in Fastq file: output/filtered/ERR1938376_filt.fastq.gz

    ## Encountered 300 unique sequences from 581 total sequences read.

    ## Not all sequences were the same length.

``` r
# Name the derep-class objects by the sample names
names(dereplicated_forward_reads) <- sample_names
```

``` r
# parameters adjusted based on recommendations for 454 data here:
# https://benjjneb.github.io/dada2/
#     faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data
dada_forward_reads <- dada(dereplicated_forward_reads,
                           err = errors_forward_reads,
                           HOMOPOLYMER_GAP_PENALTY = -1, # reduce penalty bc 454
                           BAND_SIZE = 32) # performs local alignments bc indels
```

    ## Sample 1 - 318 reads in 139 unique sequences.
    ## Sample 2 - 895 reads in 416 unique sequences.
    ## Sample 3 - 663 reads in 317 unique sequences.
    ## Sample 4 - 216 reads in 117 unique sequences.
    ## Sample 5 - 325 reads in 135 unique sequences.
    ## Sample 6 - 563 reads in 166 unique sequences.
    ## Sample 7 - 571 reads in 212 unique sequences.
    ## Sample 8 - 1267 reads in 607 unique sequences.
    ## Sample 9 - 1396 reads in 704 unique sequences.
    ## Sample 10 - 1204 reads in 697 unique sequences.
    ## Sample 11 - 1034 reads in 614 unique sequences.
    ## Sample 12 - 1079 reads in 662 unique sequences.
    ## Sample 13 - 1543 reads in 818 unique sequences.
    ## Sample 14 - 1239 reads in 744 unique sequences.
    ## Sample 15 - 1352 reads in 696 unique sequences.
    ## Sample 16 - 1162 reads in 592 unique sequences.
    ## Sample 17 - 1533 reads in 749 unique sequences.
    ## Sample 18 - 1102 reads in 528 unique sequences.
    ## Sample 19 - 956 reads in 530 unique sequences.
    ## Sample 20 - 1124 reads in 580 unique sequences.
    ## Sample 21 - 1099 reads in 523 unique sequences.
    ## Sample 22 - 1355 reads in 638 unique sequences.
    ## Sample 23 - 1263 reads in 600 unique sequences.
    ## Sample 24 - 1261 reads in 536 unique sequences.
    ## Sample 25 - 1308 reads in 659 unique sequences.
    ## Sample 26 - 1106 reads in 492 unique sequences.
    ## Sample 27 - 474 reads in 251 unique sequences.
    ## Sample 28 - 1240 reads in 551 unique sequences.
    ## Sample 29 - 970 reads in 440 unique sequences.
    ## Sample 30 - 1215 reads in 579 unique sequences.
    ## Sample 31 - 691 reads in 401 unique sequences.
    ## Sample 32 - 820 reads in 486 unique sequences.
    ## Sample 33 - 1160 reads in 419 unique sequences.
    ## Sample 34 - 1059 reads in 505 unique sequences.
    ## Sample 35 - 620 reads in 331 unique sequences.
    ## Sample 36 - 1215 reads in 668 unique sequences.
    ## Sample 37 - 1450 reads in 640 unique sequences.
    ## Sample 38 - 1219 reads in 755 unique sequences.
    ## Sample 39 - 1232 reads in 531 unique sequences.
    ## Sample 40 - 943 reads in 515 unique sequences.
    ## Sample 41 - 1392 reads in 590 unique sequences.
    ## Sample 42 - 1065 reads in 504 unique sequences.
    ## Sample 43 - 1230 reads in 633 unique sequences.
    ## Sample 44 - 1049 reads in 554 unique sequences.
    ## Sample 45 - 502 reads in 296 unique sequences.
    ## Sample 46 - 1384 reads in 667 unique sequences.
    ## Sample 47 - 1112 reads in 625 unique sequences.
    ## Sample 48 - 674 reads in 227 unique sequences.
    ## Sample 49 - 1125 reads in 291 unique sequences.
    ## Sample 50 - 1023 reads in 238 unique sequences.
    ## Sample 51 - 1189 reads in 309 unique sequences.
    ## Sample 52 - 1273 reads in 385 unique sequences.
    ## Sample 53 - 809 reads in 237 unique sequences.
    ## Sample 54 - 1098 reads in 342 unique sequences.
    ## Sample 55 - 1234 reads in 453 unique sequences.
    ## Sample 56 - 529 reads in 185 unique sequences.
    ## Sample 57 - 241 reads in 88 unique sequences.
    ## Sample 58 - 1328 reads in 440 unique sequences.
    ## Sample 59 - 723 reads in 292 unique sequences.
    ## Sample 60 - 1947 reads in 543 unique sequences.
    ## Sample 61 - 1416 reads in 527 unique sequences.
    ## Sample 62 - 1127 reads in 363 unique sequences.
    ## Sample 63 - 752 reads in 225 unique sequences.
    ## Sample 64 - 1320 reads in 530 unique sequences.
    ## Sample 65 - 1685 reads in 618 unique sequences.
    ## Sample 66 - 753 reads in 251 unique sequences.
    ## Sample 67 - 1022 reads in 362 unique sequences.
    ## Sample 68 - 1245 reads in 277 unique sequences.
    ## Sample 69 - 1531 reads in 598 unique sequences.
    ## Sample 70 - 1160 reads in 426 unique sequences.
    ## Sample 71 - 1143 reads in 395 unique sequences.
    ## Sample 72 - 1211 reads in 369 unique sequences.
    ## Sample 73 - 831 reads in 249 unique sequences.
    ## Sample 74 - 1489 reads in 470 unique sequences.
    ## Sample 75 - 1409 reads in 374 unique sequences.
    ## Sample 76 - 1075 reads in 280 unique sequences.
    ## Sample 77 - 1332 reads in 332 unique sequences.
    ## Sample 78 - 1127 reads in 360 unique sequences.
    ## Sample 79 - 1258 reads in 386 unique sequences.
    ## Sample 80 - 499 reads in 165 unique sequences.
    ## Sample 81 - 1184 reads in 256 unique sequences.
    ## Sample 82 - 1017 reads in 310 unique sequences.
    ## Sample 83 - 1109 reads in 273 unique sequences.
    ## Sample 84 - 1285 reads in 317 unique sequences.
    ## Sample 85 - 1067 reads in 294 unique sequences.
    ## Sample 86 - 1049 reads in 348 unique sequences.
    ## Sample 87 - 1056 reads in 289 unique sequences.
    ## Sample 88 - 1493 reads in 358 unique sequences.
    ## Sample 89 - 1348 reads in 336 unique sequences.
    ## Sample 90 - 742 reads in 238 unique sequences.
    ## Sample 91 - 2061 reads in 482 unique sequences.
    ## Sample 92 - 452 reads in 142 unique sequences.
    ## Sample 93 - 1499 reads in 397 unique sequences.
    ## Sample 94 - 2281 reads in 524 unique sequences.
    ## Sample 95 - 1338 reads in 389 unique sequences.
    ## Sample 96 - 1090 reads in 277 unique sequences.
    ## Sample 97 - 902 reads in 239 unique sequences.
    ## Sample 98 - 1584 reads in 489 unique sequences.
    ## Sample 99 - 1052 reads in 328 unique sequences.
    ## Sample 100 - 1277 reads in 322 unique sequences.
    ## Sample 101 - 609 reads in 130 unique sequences.
    ## Sample 102 - 2013 reads in 509 unique sequences.
    ## Sample 103 - 1252 reads in 353 unique sequences.
    ## Sample 104 - 819 reads in 187 unique sequences.
    ## Sample 105 - 277 reads in 92 unique sequences.
    ## Sample 106 - 1077 reads in 327 unique sequences.
    ## Sample 107 - 1745 reads in 577 unique sequences.
    ## Sample 108 - 1528 reads in 378 unique sequences.
    ## Sample 109 - 1191 reads in 260 unique sequences.
    ## Sample 110 - 996 reads in 277 unique sequences.
    ## Sample 111 - 841 reads in 243 unique sequences.
    ## Sample 112 - 982 reads in 377 unique sequences.
    ## Sample 113 - 462 reads in 272 unique sequences.
    ## Sample 114 - 718 reads in 267 unique sequences.
    ## Sample 115 - 581 reads in 300 unique sequences.

``` r
# check dada results
# dada_forward_reads
```

``` r
# produce the 'site by species matrix'
sequence_table <- makeSequenceTable(dada_forward_reads)
```

    ## The sequences being tabled vary in length.

``` r
# Check for and remove chimeras
sequence_table_nochim <- removeBimeraDenovo(sequence_table,
                                            method = "consensus",
                                            multithread = FALSE,
                                            verbose = TRUE)
```

    ## Identified 6 bimeras out of 1539 input sequences.

``` r
# What percent of our reads are non-chimeric?
non_chimeric_reads <- round(sum(sequence_table_nochim) / sum(sequence_table),
                            digits = 4) * 100
```

After removing chimeras, we were left with 99.97% of our cleaned reads.

``` r
# we want to export the cleaned, trimmed, filtered, denoised sequence variants
# so that we can build a phylogeny - we'll build the phylogeny outside of R
# but we need the fasta file to do so. We keep the names of each sequence as the
# sequence itself (which is rather confusing), because that's how DADA2 labels
# it's columns (e.g. 'species')
# function taken from https://github.com/benjjneb/dada2/issues/88
export_taxa_table_and_seqs <- function(sequence_table_nochim,
                                       file_seqtab,
                                       file_seqs) {
  seqtab_t <- as.data.frame(t(sequence_table_nochim)) # transpose to data frame
  seqs <- row.names(seqtab_t) # extract rownames
  row.names(seqtab_t) <- seqs # set rownames to sequences
  outlist <- list(data_loaded = seqtab_t)
  mctoolsr::export_taxa_table(outlist, file_seqtab) # write out an OTU table
  seqs <- as.list(seqs)
  seqinr::write.fasta(seqs, row.names(seqtab_t), file_seqs) # write out fasta
 }

# run the function, with the names of the files we want it to create
# and where to put them
export_taxa_table_and_seqs(sequence_table_nochim,
                           "output/sequence_variants_table.txt",
                           "output/sequence_variants_seqs.fa")
```

``` r
# assigns taxonomy to each sequence variant based on a supplied training set
# made up of known sequences
 taxa <- assignTaxonomy(sequence_table_nochim,
                       "data/training/rdp_train_set_16.fa.gz",
                       multithread = TRUE,
                       tryRC = TRUE) # also check with seq reverse compliments
```

``` r
# read in the phylogeny, which was created from the fasta exported above
# in Geneious by aligning the sequences with MAFFT and then building a
# Maximum-Likelihood tree with RAxML
 tree_in <- read_tree("output/sequence_variants_MAFFT_FastTree.newick")

# Next we want to read in the metadata file so we can add that in too
# This is not a csv file, so we have to use a slightly different syntax
# here the `sep = "\t"` tells the function that the data are tab-delimited
# and the `stringsAsFactors = FALSE` tells it not to assume that things are
# categorical variables
metadata_in <- read.table(paste0("data/metadata/",
                    "fierer_hand_bacteria_SRA_study_ERP022626_SraRunTable.txt"),
                          sep = "\t",
                          header = TRUE,
                          stringsAsFactors = FALSE,
                          row.names = 6) # sets sample IDs to row names

# Construct phyloseq object (straightforward from dada2 outputs)
phyloseq_obj <- phyloseq(otu_table(sequence_table_nochim,
                                   taxa_are_rows = FALSE), # sample-spp matrix
                         sample_data(metadata_in), # metadata for each sample
                         tax_table(taxa), # taxonomy for each sequence variant
                         phy_tree(tree_in)) # phylogeny from sequence variants
```

``` r
# bar plot of sample source (keyboard key or finger) and Abundance of sequences
plot_bar(phyloseq_obj, x = "sample_source",
         title = "Sequence abundance for all sample sources")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/sample-source-abundance-1.png) **Figure 1**: Abundance of sequences within the sapmle\_sourse(s). We see how the abundance of sequences on the each key are not all even. This sparked interest on the difference of abundance between the U and Q keys. From here a focus on the hand or finger used for each key will be compared.

``` r
# bar plot of samples from Q and U keys and Abundance of sequences
phyloseq_obj_qu <- subset_samples(phyloseq_obj,
                                  sample_source == "Ukey" |
                                    sample_source == "Qkey")
plot_bar(phyloseq_obj_qu, x = "sample_source",
         title = "Sequence abundance on Q and U keys")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-abundance-for-q-and-u-keys-1.png) **Figure 2**: Abundance of sequences within the sapmle\_sourse of the Ukey and Qkey. Since in the English language there are very few words and no common words that use the letter "Q" without "U" it was unexpected to see how the ideally more isolated Qkey having double the amount of sequence abundance. From this it will be invesitgated to see if left handed keys carry a higher abundance than right keys or the type of finger used is the cause.

``` r
# bar plot of samples from L and S keys and Abundance of sequences
phyloseq_obj_sl <- subset_samples(phyloseq_obj,
                                  sample_source == "Skey" |
                                    sample_source == "Lkey")
plot_bar(phyloseq_obj_sl, x = "sample_source",
         title = "Sequence abundance on L and S keys")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-abundance-for-s-and-l-keys-1.png) **Figure 3**: Abundance of sequences within the sapmle\_sourse of the Skey and Lkey. To separate handedness and to keep the fingers consistent, letters S and L were selected. It shows that the Lkey, which is on the right side of the keyboard, has a higher abundance. The data for the keys touched by the ringerfinger shows how either hand can led to a higher abundance. It is important to note that the Qkey had a higher abundance eventhough it is likely to be exposed or touched by the same finger that touches the Skey.

``` r
# bar plot of samples from G and H keys and Abundance of sequences
phyloseq_obj_gh <- subset_samples(phyloseq_obj,
                                  sample_source == "Gkey" |
                                    sample_source == "Hkey")
plot_bar(phyloseq_obj_gh, x = "sample_source",
         title = "Sequence abundance on G and H keys")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-abundance-for-g-and-h-keys-1.png) **Figure 4**: Abundance of sequences within the sapmle\_sourse of the Gkey and Hkey. As shown, the Gkey as on the left side has a higher abundance. The index fingers also show a higher overall abundance than the previous graphs.

``` r
# bar plot of samples from E and O keys and Abundance of sequences
phyloseq_obj_eo <- subset_samples(phyloseq_obj,
                                  sample_source == "Ekey" |
                                    sample_source == "Okey")
plot_bar(phyloseq_obj_eo, x = "sample_source",
                                  title = "Sequence abundance on E and O keys")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-abundance-for-e-and-o-keys-1.png) **Figure 5**: Abundance of sequences within the sapmle\_sourse of the Ekey and Okey. The Ekey and the Okey were selected for they are the most frequently used key on their respective side. Both the use and the abundance of the Ekey are greater than the Okey. Only one of our graphs show that a key on the right side had a greater abundance.

``` r
# box plot of left and right handed keys sequence abundance
phyloseq_table <- psmelt(phyloseq_obj)

phyloseq_table_AbundanceByKey <- phyloseq_table %>%
group_by(sample_source) %>%
summarize(sumAbundance = sum(Abundance))

phyloseq_table_AbundanceByKey$hand <- ifelse(phyloseq_table_AbundanceByKey$sample_source %in%
       c("Qkey","Wkey","Ekey","Rkey","Akey","Skey","Dkey","Fkey",
         "Gkey","Zkey","Xkey","Ckey","Vkey","Tkey","Bkey"),"L",
       ifelse(phyloseq_table_AbundanceByKey$sample_source %in%
              c("Ykey","Ukey","Ikey","Okey","Pkey","Hkey","Jkey",
                "Lkey","Nkey","Mkey"), "R", 
              ifelse(phyloseq_table_AbundanceByKey$sample_source 
                     == "Space_bar", "B","NA")
              )
       )

phyloseq_table_AbundanceByKey_LR <- phyloseq_table_AbundanceByKey %>%
  filter(hand == "L" | hand == "R")
ggplot(phyloseq_table_AbundanceByKey_LR, aes(y = sumAbundance, x = hand)) +
geom_boxplot()+ ggtitle("Comparing abundance of sequences per hand")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-abundance-between-left-and-right-hand-keys-1.png) **Figure 6**: Comparing abundance of sequences of the left and right hands. Instead of comparing letter by letter, this box plot shows the general trend that the left keys had a higher sequence abundance.

``` r
# bar plot of mean abundance for left and right handed keys

phyloseq_table_AbundanceByKey_LRB <- phyloseq_table_AbundanceByKey %>%
  filter(hand == "L" | hand == "R" | hand == "B") %>%
  group_by(hand) %>%
summarize(meanSumAbundance = mean(sumAbundance))
ggplot(phyloseq_table_AbundanceByKey_LRB, aes(x = hand,
                                              y = meanSumAbundance)) +
  geom_col() +
 ggtitle("Comparing mean abundance of sequences left and right handed keys")
```

![](Analysis_Report_01_amplicons_files/figure-markdown_github/comparing-mean-abundance-between-handed-keys-1.png) **Figure 7**: Comparing abundance of sequences for left and right handed keys and the space bar. The space bar has the highest mean abundance while the trend of the left hand having a higher sequence abundance can still be observed. The difference between the mean abundance of either hand is dwarfed when taking into condideration the keys that are touched by both. \# Discussion

Add around 2-3 pages interpreting your results and considering future directions one might take in analyzing these data.

Sources Cited
=============

Callahan,B.J. *et al.* (2016) DADA2: High-resolution sample inference from illumina amplicon data. *Nature Methods*, **13**, 581–583.

Cole,J.R. *et al.* (2008) The ribosomal database project: Improved alignments and new tools for rRNA analysis. *Nucleic acids research*, **37**, D141–D145.

McMurdie,P.J. and Holmes,S. (2013) Phyloseq: An r package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE*, **8**, e61217.
