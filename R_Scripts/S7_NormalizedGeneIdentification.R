
#' ---
#' title: "Normalize Library counts"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: allSamplesDataTable.RDS -- a table of all the samples with normalized read counts 
#'          This workflow normalizes read counts between samples to compensate for variable read depth.  
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

# load libraries
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

#' Generate load list and grouping names
#' ============================
strt <- Sys.time()
# get a list off all the files in the output folder
in.names.all <- list.files("2_output", pattern="*.rds", full.names=TRUE)
# read the load list
load.list <- read.table("input/loadlist.txt", header = FALSE, skip = 0, sep="\t",
                        stringsAsFactors = FALSE, fill=TRUE)
# set the column names
colnames(load.list) <- c("Name", "BaseName","GroupName")
# add 3 columns to the load list
load.list <- rbind(load.list,c("completeLibraryRanges","","DNA_pscAAVlib"))
# extract all the names that are in the load list
in.names.all <- in.names.all[unlist(sapply(load.list$Name, grep, in.names.all))]
# crate a datafram that holds the samples form the in.names.all and their corresponding group name based on the load list
grouping <- data.frame(Sample=gsub("-","_",gsub("found.|(2_output/)|(.rds)", "", in.names.all)),
                       Group=load.list[match(gsub("-","_",gsub("found.|(2_output/)|(.rds)", "", in.names.all)),load.list$Name),"GroupName"],
                       stringsAsFactors = FALSE)

#' Load the desired alignment files and annotating group
#' ============================

# function to load the RDS files
loadRDS <- function(in.name) {
  # read the RDS file
  this.sample <- readRDS(in.name)
  # get the name of the sample
  this.name <- gsub("-","_",gsub("found.|(2_output/)|(.rds)", "", in.name))
  # get the group of the sample
  this.group <- grouping[match(this.name,grouping$Sample),"Group"]
  # combine the sample with a data frame of the sample and the group
  mcols(this.sample) <- cbind(mcols(this.sample),
                              data.frame(Sample = this.name, Group=this.group,
                                         stringsAsFactors = FALSE))
  return(this.sample)
}

# Load all RDS files specified in 'in.names.all' into a list of R objects
all.samples <- lapply(in.names.all, loadRDS)
# Combine the list of R objects into a single GAlignmentsList object. Using unlist() to flatten the list of objects before combining
all.samples <- GAlignmentsList(unlist(all.samples))

# rename the samples into unique and R valid names
names(all.samples) <- make.names(names(all.samples), unique=TRUE)
# create a data.table with the seqnames and seq length and set the key to seqnames  
length.Table <- data.table(seqnames=names(seqlengths(all.samples)), seqlength=seqlengths(all.samples), key="seqnames")
# create a data.table from all samples and set the key to seqnames    
all.samples <- data.table(as.data.frame(all.samples), key="seqnames")
# remove unnecessary columns
all.samples[,c("strand","qwidth","cigar","njunc","end"):=NULL]
# merge all samples with the length table to get the sequence length of every sequence
all.samples <- all.samples[length.Table]
# create new columns and assign the values based on the seqnames that was splid into its parts
all.samples[, c("Category", "Protein", "Origin", 
                "Extra", "Number","GeneName") := tstrsplit(seqnames, ",", fixed=TRUE)]
# remove unnecessary columns
all.samples[, c("seqnames","Protein", "Origin", 
                "Extra", "Number") := NULL]
# cange the GeneName column to hold the actual gene name
all.samples[, GeneName := gsub("/|_","-",GeneName)]


#' Normalizing read counts to correct for variable read depth
#' ============================
# setting the key to the group
setkey(all.samples,Group)
# get rid of single read samples
all.samples <- all.samples[RNAcount>1,]
# store the read counts grouped by the group
readCounts <- all.samples[,list(GroupCount=sum(RNAcount)), by="Group"]
# normalize the group reads by the max Group reads set values between 0 and 1 where 1 is the max value
readCounts[,GroupCount:=GroupCount/max(GroupCount)]

# extract only the definitive Barcodes which have at least 2 reads per pair of barcode & LUTnr
all.samples <- all.samples[Mode=="Def"]

# set the key to Group
setkey(readCounts,Group)
# add the normalized group count to the all samples data table
all.samples <- all.samples[readCounts]
# normalize the read counts by the Group Count to correct for variable read depth
all.samples[,RNAcount:=RNAcount/GroupCount]

# set the key to the Group
setkey(all.samples,Group)
# get all the AAV samples by removing the DNA samples
total.AAV.samples <- all.samples[Group!="DNA_pscAAVlib" & Group!="DNA_pscAAVlib_Prep2" & Group!="DNA_AAVlib_DNAse_3cpc" & Group!="DNA_AAVlib_DNAse_30cpc"]
# print the number of AAV samples
print(nrow(total.AAV.samples))
# get the transported 30 cpc AAV samples
transported.AAV.samples.30cpc <- total.AAV.samples[grepl("mRNA_30cpc_SN|mRNA_30cpc_Th|mRNA_30cpc_Ctx",total.AAV.samples$Group)]
# print the number of 30 cpc AAV samples
print(nrow(transported.AAV.samples.30cpc))
# get the transported 3 cpc AAV samples
transported.AAV.samples.3cpc <- total.AAV.samples[grepl("mRNA_3cpc_SN|mRNA_3cpc_Th|mRNA_3cpc_Ctx",total.AAV.samples$Group)]
# print the number of 3 cpc AAV samples
print(nrow(transported.AAV.samples.3cpc))
# set the group name to mRNA_30cpc_Trsp
transported.AAV.samples.30cpc[,Group := "mRNA_30cpc_Trsp"]
# set the group name to mRNA_3cpc_Trsp
transported.AAV.samples.3cpc[,Group := "mRNA_3cpc_Trsp"]
# set the group name to mRNA_All
total.AAV.samples[,Group := "mRNA_All"]

# recombine all samples with their new group names
all.samples <- rbind(all.samples,total.AAV.samples,transported.AAV.samples.30cpc,transported.AAV.samples.3cpc)

# remove the data.tables that were create before
rm(total.AAV.samples,transported.AAV.samples.30cpc,transported.AAV.samples.3cpc)

# set the key of all samples to a vector of the columns
setkeyv(all.samples,c("Group","Category","GeneName","structure","start","width","Sequence","seqlength"))

# now combine all information about identical fragments that are found with different barcodes
all.samples <- all.samples[,j=list(bitScore=sum(bitScore*tCount)/sum(tCount),
                                   mismatches=median(mismatches),
                                   mCount=sum(mCount),
                                   tCount=sum(tCount),
                                   BC=paste(unique(BC), collapse = ","),
                                   Animals=paste(unique(Sample), collapse = ","),
                                   LUTnrs=paste(unique(LUTnr), collapse = ","),
                                   RNAcount=sum(RNAcount),
                                   NormCount=log2(sum(RNAcount)+1)*.N),
                           by=c("Group","Category","GeneName","structure","start","width","Sequence","seqlength")]

# adjust the start,width and seqlength from DNA based to AA based
all.samples[,start:=floor((start+2)/3)]
all.samples[,width:=ceiling((width)/3)]
all.samples[,seqlength:=ceiling(seqlength/3)]
# calculate the absolute AA position of the middle of the fragment 
all.samples[,AA:=floor(start+(width/2))]
# calculate the relative AA position of the middle of the fragment in the whole sequence
all.samples[,AAproc:=AA/seqlength*100]

#' Remove overhangs on the sequence based on the Structure annotation
#' ============================
# remove the overhangs based on the structure
all.samples[structure == "14aa", "Sequence" := substr(Sequence,3,44)]
all.samples[structure == "22aa", "Sequence" := substr(Sequence,3,68)]
all.samples[structure == "14aaG4S", "Sequence" := substr(Sequence,15,56)]
all.samples[structure == "14aaA5", "Sequence" := substr(Sequence,15,56)]

#Change the default behavior to induce start codons and Methionine
GENETIC_CODE_ALT <- GENETIC_CODE
attr(GENETIC_CODE_ALT, "alt_init_codons") <- c("TAA","TAG")

# Define the translation function
translate_function <- function(sequence) {
  as.character(Biostrings::translate(DNAString(sequence), genetic.code = GENETIC_CODE_ALT, if.fuzzy.codon = "solve"))}

# Set up parallel processing
plan(multisession, workers = detectCores())  # Adjust number of workers based on your CPU

# Apply translation function in parallel
all.samples[, Peptide := future_lapply(Sequence, translate_function)]

# change the Peptide column to character
all.samples[, Peptide := as.character(Peptide)]

# Optionally, reset the parallel plan
plan(sequential)

# print the number of samples in the all samples data table
print(nrow(all.samples))

# save the whole data table as all Samples DATA table
saveRDS(all.samples, file="0_data/allSamplesDataTable.RDS")

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()

