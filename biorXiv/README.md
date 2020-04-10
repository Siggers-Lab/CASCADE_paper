# CASCADE_paper/biorXiv
## Introduction
This folder contains an Rscript and all of the necessary data and associated annotations to recreate plots used in the main and supplementary figures of the biorXiv version of the CASCADE paper.

## Dependencies
The CASCADE pipeline and associated plotting functions were initially developed in R version 3.5.1 but should work in principle with later versions. The following R packages are required so be sure that they are installed:
| R package      | Description  |
|----------------|--------------|
|  ggplot2       | plotting package needed for the scatterplots and recruitment logos |
|  ggseqlogo     | extension of ggplot2 that actually plots the sequence logos |
|  RColorBrewer  | needed for the color palettes used in the figures |
|  cowplot       | package used to arrange figures in grids |
|  stringr       | used to facilitate string processing at certain steps in the pipeline |

## File descriptions
The following is a short description of each of the files included in this directory and how they are used by the Rscript to reconstitute the plots used in the figures from the CASCADE paper
| File                                    | Description  |
|-----------------------------------------|--------------|
|  CASCADE_figures.R                      | the R script that uses the data files and annotations below to reconstitute the CASCADE plots |
|  CXCL10_TF_sites.txt                    | positions of previously annotated TF binding sites in the CXCL10 promoter segment assayed |
|  SNP_CASCADE_Fig4.txt                   | subset of SNPs used in the Figure 4 motif grid |
|  SNP_CASCADE_exp.txt                    | experimental metadata for the SNP CASCADE array experiment |
|  SNP_SCREEN_REPRODUCIBLE_SNPs.dat       | reproducibility data for SNPs included in both the high-throughput screen and subsequenct follow-up SNP CASCADE |
|  SNP_screen_COF_exp.txt                 | experimental metadata for the differential COF recruitment screen experiments |
|  suppl_data_1_CXCL10_CASCADE.dat        | CXCL10 supplementary data from the biorXiv paper |
|  suppl_data_2_SNP_screen_COF_stats.dat  | Differential COF recruitment statistics from the biorXiv paper |
|  suppl_data_3_SNP_CASCADE.dat           | SNP CASCADE supplementary data from the biorXiv paper |

## Usage
To run the Rscript and generate the plots used in the CASCADE figures, simply run the included R script. In a UNIX-based computing environment, this may consist of loading R as a module and then calling the Rscript command as follows:
```
module load R/3.5.1
Rscript CASCADE_figures.R
```
NOTE: make sure you are running the Rscript command from the directory that contains the script itself and all of the associated data/annotations. The script does not require any arguments but expects the files described above to be located in the current working directory.

## Results
Once the script has finished running, there should be 3 additional directories each containing data and plots related to the indicated figure:

![New directories generated](screenshots/new_directories.png)

You're now ready to explore the data and plots used for the CASCADE paper!
