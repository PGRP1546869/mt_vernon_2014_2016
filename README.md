# Mt Vernon 2014/2016
Data and code for analyses contained in "Rootstock effects on scion phenotypes in a ‘Chambourcin’ experimental vineyard"

## Description 

### Morphometrics
Average morphometric data for each vine (283 vines in 2014 and 286 vines in 2016) is available in the 'mt_vernon_morpho_imagej_all_info.txt' file. R code for analyses is in the script `mt_vernon_code.R` (lines 1-137). 

### Persistent Homology
Average persistent homology data for each vine (283 vines in 2014 and 286 vines in 2016) is available in the 'mt_vernon_morpho_imagej_all_info.txt' file. R code for analyses is in the script `mt_vernon_code.R` (lines 140-292). 

### Ionomics
Ion concentrations from three leaves (young, middle, and old from a single shoot) from 288 vines for a total of 794 leaves in 2014 and 846 leaves in 2016 following outlier removal, are available in the file `mt_vernon_ionomics_2014_2016_all_data_with_vineyard_info.txt`. R code for analyses is in the script `mt_vernon_code.R` (lines 295-525).

### RNA-seq
Gene expression was examined using 28 vines with no irrigation treatment in 2016 and the normalized gene expression data is available in  `mt_vernon_rna_counts.txt`. The R package MaSigPro was used to analyse the impact of rootstock on gene expression while accounting for sampling time, since vines were sampled across 4 blocks. Experimental design information is included in the file `mt_vernon_rna_design.txt`. R code for the analyses is in the script `mt_vernon_code.R` (lines 527-574).
