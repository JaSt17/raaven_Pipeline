
#' ---
#' title: "Generate a complete library range object"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: completeLibraryRanges.rds -- a table of the AAV plasmid library so we see all fragments and their corresponding barcodes (since there are multiple barcodes for a fragment)
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

#load libraries
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))


#' Generate library rage object
#' ============================
# load the data
load("0_data/multipleContfragmentsComplete.rda") # "output.Table" output from S4: Table of all fragments reads with their corresponding barcodes
load("0_data/alignedLibraries.rda") # "allFragments.ranges" output from S3: GAlignment object with information about the fragments
load("0_data/LUTdna.rda") # "LUT.dna" LUT table with Sequence LUTnr and Structure

# change keys to make the data.table mergeable
setkey(output.Table,LUTnr)
setkey(LUT.dna,LUTnr)
# merge the tables
output.Table <- output.Table[LUT.dna,nomatch=0]
# remove unnecessary columns
output.Table[,c("Names","i.Structure"):=NULL]
# set the names of the columns
setnames(output.Table,"Sequence","fragment")
# set the key to the fragment
setkey(output.Table,fragment)
# create a new data.table with a idxFrag column which is the index of the fragments when there are multiple fragments due to multiple samples taken
range.idx <- data.table(fragment=mcols(allFragments.ranges)$Sequence, 
                        idxFrag=1:length(allFragments.ranges), key="fragment")
# only keep the fragments that are in the range.idx                         
output.Table <- output.Table[range.idx, nomatch=0, allow.cartesian=TRUE]

# create a new data.table with the fragment that are found
foundFragments.ranges <- allFragments.ranges[output.Table$idxFrag]
# remove unnecessary columns
output.Table[,c("Reads","fragment","idxFrag","Structure","LUTnr"):=NULL]
# create a new column RNAcount with the count of the fragment
output.Table[,RNAcount:=tCount]

# add the output.Table to the foundFragments.ranges
mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges), output.Table)

# save the foundFragments.ranges
saveRDS(foundFragments.ranges, file="2_output/completeLibraryRanges.rds")

devtools::session_info()