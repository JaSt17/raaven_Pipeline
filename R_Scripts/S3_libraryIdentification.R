
#' ---
#' title: "Reverse mapping of CustumArray oligos to original proteins"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: alignedLibraries.rda -- a table of the fragments with their origin protein and the LUTnr, strucutre and sequence.
#'         This workflow identifies correct fragments from the Cre-recombined AAV plasmid library and aligns them to the CustumArray ordered nucleotide fragments using Blastn.
#'         Consistant mutations in each fragment/barcode combination are also registered as is the putity of each barcode.
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(GeneGA))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(parallel))

# start timer
strt<-Sys.time()

#Override the GeneCodon function with local version containing human codons
source(file.path("functions", "GeneCodon.R"))
unlockBinding("GeneCodon", as.environment("package:GeneGA"))
assign("GeneCodon", GeneCodon, as.environment("package:GeneGA"))

#'Load sequences
#'===================
LUT.dna <- read.table("0_data/SortedFragments_all.txt",
                      header = TRUE, skip = 0, sep="\t",stringsAsFactors = FALSE, fill=TRUE)
LUT.dna <- data.table(LUT.dna)

#'Remove constitutive backbone sequences
#'===================
invisible(LUT.dna[,Sequence:=gsub("aacctccagagaggcaacg","",Sequence)])
invisible(LUT.dna[,Sequence:=gsub("cagacaagcagctaccgca","",Sequence)])
# set the sequence to upper case
invisible(LUT.dna[,Sequence:=toupper(Sequence)])
# set the key to the sequence
setkey(LUT.dna, "Sequence")
# remove duplicates
LUT.dna <- unique(LUT.dna)
# create a unqie number for each sequence
LUT.dna$LUTnr <- make.names(seq(nrow(LUT.dna)), unique=TRUE)

#'Split sequences based on linker and length 
#'===================
# 14aaG4S: GAGGCGGAGGAAGT
LUT.14aaG4S <- LUT.dna[substr(LUT.dna$Sequence,1,14) == "GAGGCGGAGGAAGT"]
# print the number of reads that match the 14aaG4S
print(nrow(LUT.14aaG4S))

# remaining sequences without 14aaG4S
LUT.remaining <- LUT.dna[!(substr(LUT.dna$Sequence,1,14) == "GAGGCGGAGGAAGT")]
# 14aaA5: CTGCTGCAGCAGCC
LUT.14aaA5 <- LUT.remaining[substr(LUT.remaining$Sequence,1,14) == "CTGCTGCAGCAGCC"]
# print the number of reads that match the 14aaA5
print(nrow(LUT.14aaA5))

# remaining sequences without 14aaA5
LUT.remaining <- LUT.remaining[!(substr(LUT.remaining$Sequence,1,14) == "CTGCTGCAGCAGCC")]
# 22aa: has 70 nucleotides and starts with CT
LUT.22aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 70L & 
                            substr(LUT.remaining$Sequence,1,2) == "CT"]
# print the number of reads that match the 22aa
print(nrow(LUT.22aa))

# 14aa: has 46 nucleotides and starts with CT
LUT.14aa <- LUT.remaining[nchar(LUT.remaining$Sequence) == 46L & 
                            substr(LUT.remaining$Sequence,1,2) == "CT"]
# print the number of sequences that match the 14aa
print(nrow(LUT.14aa))

# remove the remaining sequences because they dont match any of our 4 structures
rm(LUT.remaining)

# set the structure of the sequences
LUT.dna[LUT.dna$Sequence %in% LUT.14aaG4S$Sequence,"Structure"] <- "14aaG4S"
LUT.dna[LUT.dna$Sequence %in% LUT.14aaA5$Sequence,"Structure"] <- "14aaA5"
LUT.dna[LUT.dna$Sequence %in% LUT.22aa$Sequence,"Structure"] <- "22aa"
LUT.dna[LUT.dna$Sequence %in% LUT.14aa$Sequence,"Structure"] <- "14aa"

# save the LUT.dna with the structure the sequences and the LUTnr
save(LUT.dna,file = "0_data/LUTdna.rda")

#'Trim sequences
#'===================
LUT.14aa$Sequence <- substr(LUT.14aa$Sequence,3,44)
LUT.14aaG4S$Sequence <- substr(LUT.14aaG4S$Sequence,15,56)
LUT.14aaA5$Sequence <- substr(LUT.14aaA5$Sequence,15,56)
LUT.22aa$Sequence <- substr(LUT.22aa$Sequence,3,68)

#'Save fasta files for Bowtie alignments
#'===================
LUT.14aa.fa <- tempfile(pattern = "LUT_14aa_", tmpdir = tempdir(), fileext = "fa")
LUT.14aa.seq = ShortRead(DNAStringSet(LUT.14aa$Sequence), BStringSet(LUT.14aa$LUTnr))
writeFasta(LUT.14aa.seq,LUT.14aa.fa)

LUT.14aaG4S.fa <- tempfile(pattern = "LUT_14aaG4s_", tmpdir = tempdir(), fileext = "fa")
LUT.14aaG4S.seq = ShortRead(DNAStringSet(LUT.14aaG4S$Sequence), BStringSet(LUT.14aaG4S$LUTnr))
writeFasta(LUT.14aaG4S.seq,LUT.14aaG4S.fa)

LUT.14aaA5.fa <- tempfile(pattern = "LUT_14aaA5_", tmpdir = tempdir(), fileext = "fa")
LUT.14aaA5.seq = ShortRead(DNAStringSet(LUT.14aaA5$Sequence), BStringSet(LUT.14aaA5$LUTnr))
writeFasta(LUT.14aaA5.seq,LUT.14aaA5.fa)

LUT.22aa.fa <- tempfile(pattern = "LUT_14aaA5_", tmpdir = tempdir(), fileext = "fa")
LUT.22aa.seq = ShortRead(DNAStringSet(LUT.22aa$Sequence), BStringSet(LUT.22aa$LUTnr))
writeFasta(LUT.22aa.seq,LUT.22aa.fa)

#'Build Bowtie index
#'===================
# read in the sequneces from the original fasta file
seqs.original <- readFasta("input/DNA-lib_RetrogradeTransport.fasta")

# translate the sequences form all RetrogradeTransport sequences to amino acids
seqs.AA <- Biostrings::translate(sread(seqs.original), genetic.code=GENETIC_CODE, 
                                 if.fuzzy.codon="error")

# load the function to convert amino acids to DNA
source("functions/AAtoDNA.R")

# convert the amino acids back to DNA but using human specific codons
seqs.optimized = ShortRead(DNAStringSet(sapply(seqs.AA, function(x) AAtoDNA(x, species="hsa"))), 
                           BStringSet(gsub("([ ])", "_", ShortRead::id(seqs.original))))

# write the optimized sequences to a fasta file
bowtie.fasta <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(seqs.optimized,bowtie.fasta)

# create a temporary file for the bowtie index
bowtie.idx <- tempfile(pattern = "IDX_bowtie_", tmpdir = tempdir(), fileext = "")

sys.out <-  system(paste("bowtie2-build",bowtie.fasta,bowtie.idx, "2>&1",  sep = " "), 
                   intern = TRUE, ignore.stdout = FALSE) 


#' Align fragments to reference
#' ============================
# Function to align sequences and process the results
align_and_process <- function(fa_file, lut, bowtie_idx) {
  # Create a temporary file name for Bowtie output
  name.bowtie <- tempfile(pattern = "bowtie_", tmpdir = tempdir(), fileext = "")
  
  # Run Bowtie2 alignment
  sys.out <- system(paste("bowtie2 --non-deterministic --threads ", detectCores(),
                          " --very-sensitive -f -a",
                          " -x ", bowtie_idx, " -U ", fa_file, " -S ", 
                          name.bowtie, ".sam 2>&1", sep = ""), 
                    intern = TRUE, ignore.stdout = FALSE)
  
  # Run samtools to convert the SAM file to a BAM file
  command.args <- paste("view -@ ", detectCores(), " -Sb ", name.bowtie, ".sam > ",
                        name.bowtie, ".bam", sep = "")
  system2("samtools", args = command.args, stdout = TRUE, stderr = TRUE)
  
  # Run samtools to sort the BAM file
  command.args <- paste("sort -@ ", detectCores(), " ", name.bowtie, ".bam -o ",
                        name.bowtie, "_sort.bam", sep = "")
  system2("samtools", args = command.args, stdout = TRUE, stderr = TRUE)
  
  # Read the sorted BAM file into a GAlignments object
  frag_ranges <- readGAlignments(paste(name.bowtie, "_sort.bam", sep = ""), use.names = TRUE)
  
  # Return the results as a list
  list(
    frag_ranges = frag_ranges,
    total = length(names(frag_ranges)),
    unique = length(unique(names(frag_ranges))),
    lut_unique = length(unique(lut$Sequence))
  )
}
# run the function for the different libraries
# Align and process 14aa sequences
results_14aa <- align_and_process(LUT.14aa.fa, LUT.14aa, bowtie.idx)
# Align and process 14aaG4S sequences
results_14aaG4S <- align_and_process(LUT.14aaG4S.fa, LUT.14aaG4S, bowtie.idx)
# Align and process 14aaA5 sequences
results_14aaA5 <- align_and_process(LUT.14aaA5.fa, LUT.14aaA5, bowtie.idx)
# Align and process 22aa sequences
results_22aa <- align_and_process(LUT.22aa.fa, LUT.22aa, bowtie.idx)

# Annotate the aligned sequences with their structure
mcols(results_14aa$frag_ranges)$structure <- "14aa"
mcols(results_22aa$frag_ranges)$structure <- "22aa"
mcols(results_14aaA5$frag_ranges)$structure <- "14aaA5"
mcols(results_14aaG4S$frag_ranges)$structure <- "14aaG4S"

# Merge all the aligned sequences into one object
allFragments.ranges <- append(results_14aa$frag_ranges, results_22aa$frag_ranges)
allFragments.ranges <- append(allFragments.ranges, results_14aaA5$frag_ranges)
allFragments.ranges <- append(allFragments.ranges, results_14aaG4S$frag_ranges)

# Annotate the merged sequences with their LUT number
mcols(allFragments.ranges)$LUTnr <- names(allFragments.ranges)
# Set the key for the LUT data frame
setkey(LUT.dna, LUTnr)
# Annotate the merged sequences with their sequence from the LUT based on the LUT number
mcols(allFragments.ranges)$Sequence <- LUT.dna[mcols(allFragments.ranges)$LUTnr]$Sequence

# Save the merged and annotated sequences to a file
save(allFragments.ranges, file = "0_data/alignedLibraries.rda")

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()

