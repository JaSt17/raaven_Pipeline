
#' ---
#' title: "Library fragment alignment"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: multipleContfragmentsComplete.rda -- a table holding |BC|tCount|Mode|mCount|mismatches|bitScore|LUTnr|Structure
#'         This workflow identifies correct fragments from the Cre-recombined AAV plasmid library and aligns them to the CustumArray ordered nucleotide fragments using Blastn.
#'         Consistant mutations in each fragment/barcode combination are also registered as is the putity of each barcode.
#' ---


suppressPackageStartupMessages(library(knitr)) 

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(kableExtra))

strt1<-Sys.time()

#' Load the trimmed reads
#' ============================

# loading the LUT and the fragments
load("0_data/LUTdna.rda")

# please select the barcode and fragment files that you want to analyze (Prep1 or Prep2)
filename <- "Prep1"

# change the filename to the correct path
path <- paste0("0_data/barcode_fragment/",filename,"_DNA_pscAAVlib/")

# loading the fragments and the barcodes 
fragments.file <- paste0(path,"fragments_AAVlibrary_complete.fastq.gz")
barcodes.file <- paste0(path,"barcodes_AAVlibrary_complete.fastq.gz")

# read the fragments and the barcodes into fastq files
reads.trim <- readFastq(fragments.file)
reads.BC <- readFastq(barcodes.file)


#' Make CustomArray reference index for Blast
#' ============================
#+ Making Bowtie index.......

# convert the seq form the LUT to a fasta file for downstream analysis
LUT.fa <- tempfile(pattern = "LUT_", tmpdir = tempdir(), fileext = ".fa")
LUT.seq = ShortRead(DNAStringSet(LUT.dna$Sequence), BStringSet(1:length(LUT.dna$LUTnr)))
writeFasta(LUT.seq,LUT.fa)


#'Save unique fragments as fasta file
#'===================
unique.reads <- unique(sread(reads.trim))
# converting the fragments to a fasta file for downstream analysis
unique.reads <- ShortRead(DNAStringSet(unique.reads), BStringSet(1:length(unique.reads)))
fragments.unique.fa <- tempfile(pattern = "FragUnique_", tmpdir = tempdir(), fileext = ".fa")
writeFasta(unique.reads,fragments.unique.fa)


#'Align against the library using blast
#'===================
# create a temporary file for the blast database and the output
blast.db <- tempfile(pattern = "blastDB_", tmpdir = tempdir(), fileext = ".db")
blast.out <- tempfile(pattern = "blastOut_", tmpdir = tempdir(), fileext = ".txt")

# create the blast database using the LUT sequnces
sys.out <-  system(paste("makeblastdb -in ", LUT.fa,
                         " -out ",blast.db," -dbtype nucl -title LUT -parse_seqids 2>&1",  sep = ""),
                   intern = TRUE, ignore.stdout = FALSE) 

# print the size of the database
system(paste("blastdbcmd -db ", blast.db, " -info", sep=" "))
# print the number of unique fragments
system(paste("grep -c '>' ", fragments.unique.fa, sep=" "))

# Determine the number of cores available
num_threads <- detectCores()

# Construct the blastn command
command <- paste0(
  "blastn -max_target_seqs 25 -word_size 7",
  " -num_threads ", num_threads,
  " -outfmt 10 -db ", shQuote(blast.db),
  " -query ", shQuote(fragments.unique.fa),
  " > ", shQuote(blast.out), " 2>&1"
)
# run the blastn look for all fragments in the db created by the sequneces in the LUT
sys.out <- system(command, intern = TRUE, ignore.stdout = FALSE)

# write the output to a file
system(paste("gzip -c ", blast.out, " > ./0_data/", filename, "_blastOutput.csv.gz", sep=""))

# create a data.table from the blast output
table.blastn <- data.table(scan(file=paste("./0_data/", filename, "_blastOutput.csv.gz", sep=""), what="character", sep=";") , keep.rownames=FALSE, key="V1")

# if there are warnings in the blast output, remove them and store them in a separate data.table
if (length(grep("Warning",table.blastn$V1)) != 0) {
  warnings.out <- unique(table.blastn[grep("Warning",table.blastn$V1),])
  table.blastn <- table.blastn[-grep("Warning",table.blastn$V1),]
  setnames(warnings.out,"V1", c("blastn Warnings"))
}

# create new columns and fill them by splitting the V1 column
table.blastn[,c("Reads","Sequence","identity","alignmentLength","mismatches",
                "gapOpens", "q_start", "q_end", "s_start", "s_end", 
                "evalue","bitScore") := tstrsplit(V1,",",fixed=TRUE),]

# replacing the numeric values with the actual sequence data from unique.reads and LUT.seq
table.blastn[,Reads:= as.character(sread(unique.reads[as.integer(Reads)]))]
table.blastn[,Sequence:= as.character(sread(LUT.seq[as.integer(Sequence)]))]


#'Create full table of all fragments that mached the LUT with their blastn results
#'===================
# set the key to the sequence in both tables
setkey(table.blastn,Sequence)
setkey(LUT.dna,Sequence)

# keep only the rows that have a match in the LUT
table.blastn<- table.blastn[LUT.dna, nomatch=0]
# remove the sequence column
table.blastn[,c("V1","identity","alignmentLength","gapOpens", "q_start", 
                "q_end", "s_start", "s_end", "evalue","Sequence","Names"):=NULL]

# change the data types of the columns
table.blastn[,bitScore:= as.numeric(bitScore)]
table.blastn[,mismatches:= as.numeric(mismatches)]

# set the key to the reads and the LUTnr
setkeyv(table.blastn,c("Reads","LUTnr"))
# make sure that the fragments are only aligned once to the reference in the top ten matches
setorder(table.blastn,Reads,LUTnr,-bitScore)
# have only unique pairs of reads and LUTnr
table.blastn <- unique(table.blastn, by=c("Reads","LUTnr"))

# create a full table with the fragments and the barcodes
full.table <- data.table(Reads=as.character(sread(reads.trim)),
                         BC=as.character(sread(reads.BC)),
                         key="Reads")

# selecting only the top hit for each read so we now have only one read for each LUTnr
table.blastn.topHit <- table.blastn[table.blastn[, .I[which.max(bitScore)], by="Reads"]$V1]
# merge the full table with the top hit table so we have only the best hit reads
full.table <- full.table[table.blastn.topHit,nomatch=0]

# get the number of rows in the full table corresponding to the number of reads with barcodes
all.reads <- nrow(full.table)

# calculate the percentage of reads that have been aligned to our db 
print(paste("Alignment percentage:", percent(nrow(full.table)/all.reads)))


#' Starcode based barcode reduction
#' ============================
# create a temporary file for the output
out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")
# create the command for the starcode that clusters barcodes which are only one mismatch away from each other together
command.args <- paste("-c ",barcodes.file," | starcode -t ",detectCores()-1," --print-clusters -d",
                      1," -r5 -q -o ", out.name.BC.star, sep = "")
# run starcode 
system2("gunzip", args=command.args, stdout=TRUE, stderr=TRUE)
# read the output from starcode into a data.table
table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                                     stringsAsFactors = FALSE, fill=FALSE),keep.rownames=TRUE, key="rn")
# delete the V2 column
table.BC.sc[,V2 := NULL]
# split the V3 column into separate columns
table.BC.sc <- table.BC.sc[, strsplit(as.character(V3),",",fixed=TRUE), by=rn]
# calculate the number of barcodes that have been dropped
SC.droppedBC <- length(unique(sread(reads.BC))) - length(unique(table.BC.sc$V1))
print(paste("Dropped BCs in Starcode:", SC.droppedBC))
# rename the columns in the table to BC and scBC so we have the barcode and its corresponding starcode reduced barcode
setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))
# save the table of barcodes and their starcode reduced versions
DNA_pscAAVlib <- table.BC.sc
save(DNA_pscAAVlib, file=paste("0_data/", filename, "_scBC_DNA_pscAAVlib.rda", sep=""))


#' Replacing barcodes with Starcode reduced versions
#' ============================
# set the keys in both tables
setkey(full.table,BC)
setkey(table.BC.sc,BC)
# merge the full table with the starcode reduced table based on their BC
full.table <- full.table[table.BC.sc,nomatch=0]

# change the names of the columns
setnames(full.table,c("BC","scBC"),c("oldBC","BC"))
# change the key to the new barcode
setkey(full.table,BC)

# calculate the number of original barcodes and the number of unique barcodes after the reduction
RetainedBC <- length(unique(full.table$oldBC))
scBC <- length(unique(full.table$BC))
print(paste("Original unique barcodes:", RetainedBC))
print(paste("SC reduced unique barcodes:", scBC))

# remove the old barcode column
invisible(full.table[,oldBC:=NULL])


#' Splitting reads into single-read and multi-read barcodes
#' ============================
# order the full table by the barcode and change the mismatches to numeric
full.table <- full.table[order(full.table$BC),]
full.table[,mismatches:= as.numeric(mismatches)]

# split the full table into single and multi-read barcodes
temp.table.single <- full.table[full.table[, .I[.N == 1], by="BC"]$V1]
temp.table.multi <- full.table[full.table[, .I[.N > 1], by="BC"]$V1]

# printig the number of reads we have and the number of single reads
print("Utilized reads.......")
print(nrow(full.table))
print("Whereof single reads.......")
print(nrow(temp.table.single))

# set mcount and tcount to 1 for the single-read table
temp.table.single[,c("mCount","tCount"):=1]
# set the mode to Amb for ambiguous
temp.table.single$Mode <- "Amb"

# setting a key vector for the multi-read table
setkeyv(temp.table.multi,c("BC","LUTnr"))
# set the name of the mode to Def for definitive
temp.table.multi$Mode <- "Def"

# add three new columns to the multi-read table the mean bitScore, median mismatches and the total number of reads based on pairs of the same BC and LUTnr
temp.table.multi[,c("bitScore","mismatches" ,"tCount"):= list(mean(bitScore),
                                                              median(mismatches), .N), 
                 by=key(temp.table.multi)]
# tcount is the number of reads for each pair of barcode and LUTnr

# remove duplicates since we now have the tcount and the average bitScore and median mismatches we now have only uqniue combination of BC and LUTnr
temp.table.multi <- unique(temp.table.multi)
# now we only have unique pairs of BC and LUTnr in the multi-read table


#' Splitting multi-read barcodes into clean and chimeric
#' ============================
#+ Splitting Clean Reads.......

# set the key to the barcode
setkeyv(temp.table.multi,"BC")

# split the multi-read barcodes into once that have only one LUTnr and those that have more than one LUTnr
temp.table.multi.clean <- temp.table.multi[temp.table.multi[, .I[.N == 1], by="BC"]$V1] # clean since 1 to 1 relation between BC and LUTnr
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[.N > 1], by="BC"]$V1] # chimeric since 1 to many relation between BC and LUTnr

# set the mCount to the tCount in barcodes with only one LUTnr so mCount= number of reads for this pair of BC and LUTnr
temp.table.multi.clean[,mCount:=tCount]

# print the number of clean and chimeric multi-read barcodes
print("Clean multi-read barcodes.......")
print(nrow(temp.table.multi.clean))
print("Chimeric multi-read barcodes.......")
print(length(unique(temp.table.multi$BC)))


#' Calculate consensus alignment of chimeric barcodes
#' ============================
# set the key to the barcode
setkey(temp.table.multi,"BC")
# set the mCount to the number of reads for each pair of BC and LUTnr
temp.table.multi[, mCount:=tCount]

# now set the tCount to the sum of the tCounts of all entries with the same barcode (not considering the LUTnr) (tCount is the total number of reads for each barcode)
temp.table.multi[, "tCount":=sum(tCount), by="BC"]

# delete the columns that are not needed
temp.table.multi[,c("LUTnr","bitScore","mismatches","Structure"):=NULL]

# set the key to the Reads
setkey(temp.table.multi,"Reads")
# set the key of table.blastn to the Reads
setkey(table.blastn,"Reads")

# join the table.blastn with the temp.table.multi based on the Reads. We allow for multiple matches and only keep the ones that have a match in both tables
temp.table.multi <- temp.table.multi[table.blastn, nomatch=0, allow.cartesian=TRUE]

# set the key to a vector of the barcode and the LUTnr
setkeyv(temp.table.multi,c("BC","LUTnr"))

# change the columns so we get the max bitScore, median mismatches and the sum of the mCounts for each pair of BC and LUTnr
temp.table.multi[,c("bitScore","mismatches" ,"mCount"):= list(max(bitScore),
                                                              median(mismatches), 
                                                              sum(mCount)), by=key(temp.table.multi)]
# mCount is the number of reads for each pair of barcode and LUTnr if we not only take the best hit from blast for each read

# delete duplicates so we now only have unique combination of BC and LUTnr with the highest bitScore
temp.table.multi <- unique(temp.table.multi, by=c("BC","LUTnr"))

# set the key to the barcode
setkeyv(temp.table.multi,"BC")
# Select only BC LUTnr pairs with the maximum mCount for which means the maximum number of reads that would be possible for this barcode
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[mCount == max(mCount)], 
                                                      by=key(temp.table.multi)]$V1] 
# If there are multiple BC LUTnr pairs having the same number of reads we chose the barcode LUTnr pait with the highest bitScore from the blastn search
temp.table.multi <- temp.table.multi[temp.table.multi[, .I[which.max(bitScore)], 
                                                      by=key(temp.table.multi)]$V1] 

# if the barcode LUTnr combination has the highst bitScore and only has one mCount then it is a clean barcode
temp.table.multi[temp.table.multi$mCount==1]$Mode <- "Amb"

# create a combined table with the clean cimeric and clean barcodes
temp.table.multi.consensus <- rbind(temp.table.multi, temp.table.multi.clean)

print(paste("Total number of definitive Barcodes:", 
            length(grep("Def", temp.table.multi.consensus$Mode))))
print(paste("Total number of ambiguous Barcodes:", 
            length(grep("Amb", temp.table.multi.consensus$Mode))))
print(paste("Total number of single-read Barcodes:", 
            nrow(temp.table.single)))

# create the final output table with the single and multi-read barcodes
output.Table <- rbind(temp.table.multi.consensus,temp.table.single)
save(output.Table, file=paste("0_data/", filename,"_multipleContfragmentsComplete.rda")

print("Total analysis time:")
print(Sys.time()-strt1)
devtools::session_info()
