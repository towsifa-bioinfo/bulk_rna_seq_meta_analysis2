# Install pak package manager
install.packages("pak")

# Install packages using pak
pak::pkg_install(c("tidyverse", "DESeq2", "VennDiagram", "ggVennDiagram", "ggrepel", "biomaRt", "ensembldb", "AnnotationHub"))

# Install metavolcano package
# https://github.com/csbl-usp/MetaVolcanoR
# Install devtools if needed
install.packages("devtools")

# Install MetaVolcanoR
devtools::install_github("csbl-usp/MetaVolcanoR")

# Test installed packages
library(tidyverse)
library(DESeq2)
library(VennDiagram)
library(ggVennDiagram)
library(MetaVolcanoR)
library(biomaRt)
library(ensembldb)
library(AnnotationHub)
