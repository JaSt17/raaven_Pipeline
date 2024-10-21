
#' ---
#' title: "Extraction of Barcodes and gene fragments"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: barcode_***.fastq.gz, fragments_***.fastq.gz -- fastq.gz files with the paired barcodes and fragments from the given Library.
#'        This workflow extracts the barcodes and fragments from the sequencing files and saves them as fastq.gz files.
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(kableExtra))

config <- read.table("input/config.txt", header = FALSE, skip = 0, sep="\t",
                     stringsAsFactors = FALSE, fill=TRUE)
colnames(config) <- c("Parameter", "Value")

#'Sequencing files
#'===================
# load the data directory and the input files
dataDir <- config$Value[1]
in.name.P5 <- file.path(dataDir, config$Value[2])
in.name.P7 <- file.path(dataDir, config$Value[3])
name.out <- config$Value[4]

strt<-Sys.time()

# Counting reads and displaying the number of reads
output.Reads <- as.integer(system(paste("gunzip -c ",shQuote(gsub("([\\])", "", in.name.P5)),
                                              " | echo $((`wc -l`/4)) 2>&1", sep = ""), intern = TRUE, 
                                        ignore.stdout = FALSE))
print(paste("Utilized Reads:", output.Reads))


#' Extraction of barcodes
#' ============================

# create a temporary file for the output
out.name.P5 <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

# run the bbduk2 command
sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=18 mink=18 hammingdistance=2 findbestmatch=t ",
                        "rcomp=f findbestmatch=f qhdist=1 minavgquality=0 maxns=0 minlength=18 ",
                        "maxlength=22 threads=", detectCores()," in=", shQuote(in.name.P5), 
                        " out=", out.name.P5," lliteral=", "GGCCTAGCGGCCGCTTTACTT",
                        " rliteral=", "ATAACTTCGTATAATGTATGC",
                        " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) 

# set the output as the new input
in.name.P5 <- out.name.P5

# remove the output
rm(sys.out)

# read the barcodes
reads.BC <- readFastq(in.name.P5)
sread(reads.BC)
print(paste("Total number of found barcode reads:", length(reads.BC)))


#' Extraction of fragments
#' ============================
#+ Extracting fragments.......

# create a temporary file for the output
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
# run the bbduk2 command
command.args <- paste("overwrite=true k=18 mink=18 rcomp=f qhdist=1 maskmiddle=t",
                      " hammingdistance=2 findbestmatch=t minlength=38 maxlength=78 ordered=t ",
                      "threads=", detectCores(), " in=", in.name.P7, " out=", out.name.P7,
                      " lliteral=", "AGCAACCTCCAGAGAGGCAACG",
                      " rliteral=", "CAGACAAGCAGCTACCGCAGAT", sep = "")

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) 

# remove the output
rm(sys.out)

# set the input as the output
in.name.P7 <- out.name.P7

# read the fragments
reads.Frag <- readFastq(in.name.P7)
sread(reads.Frag)
print(paste("Total number of found fragment reads:", length(reads.Frag)))

# create a temporary file for the output of the fragments and barcodes 
out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P5_singlet <- tempfile(pattern = "P5_singlet_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7_singlet <- tempfile(pattern = "P7_singlet_", tmpdir = tempdir(), fileext = ".fastq.gz")

# run the pairfq command to pair the barcodes and fragments
command.args <- paste("makepairs -c 'gzip' -f ", in.name.P5," -r ", in.name.P7,
                      " -fp ", out.name.P5, " -rp ", out.name.P7, " -fs ",
                      out.name.P5_singlet, " -rs ", out.name.P7_singlet,
                      " --stats 2>&1", sep = "")

sys.out <- system2("pairfq", args=command.args, stdout=TRUE, stderr=TRUE)

# save the barcodes and fragments as fastq.gz files in 0_data
system(paste("mv ", out.name.P5, " ./0_data/barcodes_", name.out, ".fastq.gz", sep=""))
system(paste("mv ", out.name.P7, " ./0_data/fragments_", name.out, ".fastq.gz", sep=""))

# clean up temporary files
unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE)

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()

