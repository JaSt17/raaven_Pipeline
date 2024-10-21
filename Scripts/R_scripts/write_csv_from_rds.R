
#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

# load libraries
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

# Load the rds file into a dataframe
all.samples <- readRDS("0_data/allSamplesDataTable.RDS")

all.samples$Peptide <- as.character(all.samples$Peptide)

# Write the dataframe to a csv file
write.csv(all.samples, "0_data/allSamplesDataTable.csv", row.names = FALSE)