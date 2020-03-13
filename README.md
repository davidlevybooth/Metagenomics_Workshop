# Visualization of Metagenomic and Metatranscriptomic Data in R | ECOSCOPE

## Overview

Metagenomics and metatranscriptomics offer unparalleled insight into the taxonomic and functional diversity present across a multitude of environments. Yet the vast amounts of data these tools produce is challenging to effectively analyze and visualize.

We will analyze a unique thermal swamp metagenome and metatranscriptome dataset to develop strategies for effective visual communication.

![Thermal Swamp Metagenomics](images/metagenomics.header.png)

This workshop is designed around the application of whole-genome phylogenetic analysis to organize the diverse data types that we can generate from metagenomics and metatranscriptomics analysis. 

You will learn: 

1. Manipulating and visualizing whole-genome phylogenies with <i>ggtree</i> 
2. Mapping genome assembly statistics to phylogenetic trees with <i>ggplot2</i>
3. Processing count data from <i>Meiji</i> and mapping to phylogenetic trees
4. Genomic functional annotation and pathway analysis with <i>pathview</i>
5. Mapping metatranscriptomic count data and coverage to circular genome maps with <i>circlize</i>
6. Bonus: Mapping metatranscriptomic count data and coverage to linear contig maps

This workshop assumes some prior experience with R (such as that in our [Introduction to R workshop](https://github.com/EDUCE-UBC/workshops_R/tree/master/intro_R_2hr)). All code presented in this workshop is contained in the “metagenomics.R” R script file. 

You can follow along with the workshop by selecting the relevant line(s) of code and pressing ctrl+enter on Windows, cmd+enter on MacOS, or using the “Run” button in the upper right of your script window to execute it.

Note that the majority of this workshop uses base R functions. No experience with Tidyverse is required. 

Please call over the instructor or a TA if you have any problems with executing the code or if you need any help with the exercises. 

## Prior to the workshop

# Setup Instructions

Please come to the workshop with your laptop setup with the required software and data files as described in our setup instructions.

# Making an RStudio Project

Projects allow you to divide your work into self-contained contexts.

Let’s create a project to work in.

In the top-right corner of your RStudio window, click the “Project: (None)” button to show the projects
drop-down menu. Select “New Project…” > “New Directory” > “New Project.” Under directory name, input
“metagenomics” and choose a parent directory to contain this project on your computer.

### Installing Packages

Many R packages can be downloaded and installed from The Comprehensive R Archive Network (CRAN).


You can install each package with the following code: 

```{r eval = FALSE}
install.packages("BiocManager",
                   "devtools",
                   "RColorBrewer",
                   "ape", 
                   "tidyverse",
                   "ggplot2",
                   "phytools",
                   "stringr", 
                   "reshape2",
                   "cowplot",
                   "grid",
                   "circlize",
                   "dplyr")

```

To make things a bit more difficult, there are other source of R packages, such as Bioconductor, which provides tools for the analysis and comprehension of high-throughput genomic data.

```{r eval = FALSE}
BiocManager::install("pathview")}
```

##Install ggtree

And finally, we have to install the developmental version of _ggtree_ from [GitHub](https://github.com/YuLab-SMU/ggtree). **This is very important to ensure colored ranges are applied correctly** (the Bioconductor version does not apply colors to trees correctly)

```{r eval = FALSE}
library(devtools)
install_github("YuLab-SMU/ggtree")
```

#Loading Packages 

Next, we should check if all packages can be loaded properly. If a package does not load because of a missing package, please install the package (as shown above) and try again.  

```{r, message=FALSE, warning=FALSE, eval = FALSE}
# Suite of packages for data manipulation and visualization
library(tidyverse)
# Working with trees
library(ape)
library(phytools)
# Plotting Grammar of Graphics
library(ggplot2)
# Tree visualization
library(ggtree)
# String manipulation
library(stringr)
# Data structure manipulation
library(reshape2)
# Multi-panel figure layouts
library(cowplot)
library(grid)
# Color maniputation
library(RColorBrewer)
# Pathway analysis with KEGG Orthologs
library(pathview)
#plot circular graphs and contig maps
library(circlize)
```

#Downloading Workshop Material

The material for this workshop is hosted on [GitHub](https://github.com/davidlevybooth/Metagenomics_Workshop).

Follow the above link and click the green "Clone or download"" button to download all files as a .zip archive. Unzip the archive to access files. 

# Loading custom functions

I have written some helpful functions to replace some time-consuming workflows. Please load them from the source R script like this: 

```{r}
#Custom functions for metagenomics and metatranscriptomics
source("scripts/metagenome.utilities.R")
```

That's it! 

Please ensure all packages are installed and up to date, and that you have all data downloaded prior to the start of the workshop. If you require any assistance, please contact David Levy-Booth (dlevyboo@mail.ubc.ca). 

We will also spend a few minutes prior to the workshop ensuring that all attendees have completed setup and installation.
