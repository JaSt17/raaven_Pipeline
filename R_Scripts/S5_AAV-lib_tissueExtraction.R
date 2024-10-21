
#' ---
#' title: "Barcoded extraction and reduction from RNA samples"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output: found.***.rds -- a table of the found fragments in the RNA samples given by the list in the input/loadlist.txt. 
#'         This workflow extracts the barcode and then reducedes them using the starcode algorithm. They get 
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

#load libraries
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(beanplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(biovizBase))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))


#' Analyze tissue RNA
#' ============================
# start the timer
strt <- Sys.time()
# load the data
load("0_data/multipleContfragmentsComplete.rda") # "output.Table" output from S4: Table of all fragments reads with their corresponding barcodes
load("0_data/alignedLibraries.rda") # "allFragments.ranges" output from S3: GAlignment object with information about the fragments
load("0_data/LUTdna.rda") # "LUT.dna" LUT table with Sequence LUTnr and Structure

# load the list of the RNA samples that should be analyzed
load.list <- read.table("input/loadlist.txt", header = FALSE, skip = 0, sep="\t",
                        stringsAsFactors = FALSE, fill=TRUE)
# set the path to the fastq files
dataDir <- "0_data/fastq_files"
colnames(load.list) <- c("Name", "BaseName","GroupName")
# create a empty log table
log.table <- data.table(Name="Name",
                        Reads=NA,
                        BCs=NA,
                        allBCs=NA, # all barcodes
                        scBCs=NA,
                        SCdroppedBC=NA) # number of barcodes dropped by Starcode) # Starcode reduced barcodes

# function to analyze tissue by their index number
analyzeTissue <- function(indexNr) {
  
  # Extract the base name for the file corresponding to the provided index
  baseName <- load.list$Name[indexNr]
  

  #' If necessary, split the base name into components using '/' as a delimiter to separate subdirectories
  #' ============================
  nameComponents <- unlist(strsplit(baseName, "/"))
  # Remove any NA values from the name components
  nameComponents <- nameComponents[!is.na(nameComponents)]
  # Determine the input file paths based on the number of name components
  if (length(nameComponents) == 2) {
    # If the name has two components, construct the path to the files within the subdirectory
    in.files <- list.files(
      path = file.path(gsub("([\\])", "", dataDir), nameComponents[1]), 
      pattern = paste0(nameComponents[2], "*"), 
      full.names = TRUE
    )
  } else {
    # If the name has a single component, search for files directly within dataDir
    in.files <- list.files(
      path = gsub("([\\])", "", dataDir), 
      pattern = paste0(nameComponents[1], "*"), 
      full.names = TRUE
    )
  }
  
  # Identify the files corresponding to P5 (e.g., _1 in the name) because the barcodes are in the P5 reads
  in.name.P5 <- in.files[grep("_1", in.files)]

  # add the name to the log table and set the output name
  log.table$Name <- load.list$Name[indexNr]
  log.table$Reads <- length(readFastq(in.name.P5))
  name.out <- log.table$Name


  # Extraction of barcodes
  # ============================
  # create a temporary file for the output
  out.name.BC <- tempfile(pattern = "BC_", tmpdir = tempdir(), fileext = ".fastq.gz")

  # run the bbduk2 command to extract barcodes fromt the reads
  sys.out <- system(paste("~/bbmap/bbduk2.sh overwrite=true k=12 mink=12 hammingdistance=2 ",
                          "findbestmatch=t trd=t rcomp=f skipr2=t findbestmatch=f qhdist=0 ",
                          "minavgquality=0 ordered=t maxns=0 minlength=18 maxlength=22 threads=", 
                          detectCores()," in=", shQuote(in.name.P5), " out=", out.name.BC,
                          " lliteral=", "GGCCTAGCGGCCGCTTTACTT", " rliteral=", "ATAACTTCGTATA",
                          " 2>&1", sep = ""), intern = TRUE, ignore.stdout = FALSE) 
                          
  # save the Barcode extraction result in the log table
  log.table$BCs <- strsplit(sys.out[grep("Result:",sys.out)],split = "\t")[[1]][2]
  # read the barcodes
  reads.BC <- readFastq(out.name.BC)
  # create a table with the barcodes and the ID
  barcodeTable <- data.table(ID=as.character(ShortRead::id(reads.BC)), 
                            BC=as.character(sread(reads.BC)), key="BC")


  # Starcode based barcode reduction
  # ============================
  # create a temporary file for the output
  out.name.BC.star <- tempfile(pattern = "BCsc_", tmpdir = tempdir(), fileext = ".txt")
  # run the starcode command to reduce the barcodes
  system(paste("gunzip -c ",out.name.BC," | starcode -t ",detectCores()-1," --print-clusters -d",
              1," -r5 -q -o ", out.name.BC.star, " 2>&1", sep = ""), 
        intern = TRUE, ignore.stdout = FALSE)
  # create a table with the reduced barcodes
  table.BC.sc <- data.table(read.table(out.name.BC.star, header = FALSE, row.names = 1, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=FALSE),keep.rownames=TRUE, key="rn") 
  # remove the second column
  table.BC.sc[,V2 := NULL]
  # split the third column by comma in respect to the grouping by their Starcode reduced barcode
  table.BC.sc <- table.BC.sc[, strsplit(as.character(V3),",",fixed=TRUE), by=rn]
  # set the column names to BC and scBC
  setnames(table.BC.sc,c("V1","rn"),c("BC","scBC"))


  # Replacing barcodes with Starcode reduced versions
  # ============================
  # set the key to BC
  setkey(table.BC.sc,BC)
  # subset the barcodeTable with the reduced barcodes 
  barcodeTable <- barcodeTable[table.BC.sc,nomatch=0]
  # set the column names to oldBC and BC
  setnames(barcodeTable,c("BC","scBC"),c("oldBC","BC"))
  # set the key to BC
  setkey(barcodeTable,BC)

  # save the number of unqiue old barcodes and remaining barcodes in the log table
  log.table$allBCs <- length(unique(barcodeTable$oldBC))
  log.table$scBCs <- length(unique(barcodeTable$BC))
  log.table$SCdroppedBC <- log.table$allBCs - log.table$scBCs

  # delete the old barcodes column
  invisible(barcodeTable[,oldBC:=NULL])
  # set the key of the output table to BC
  setkey(output.Table,"BC")

  # reorganize the barcodeTable to to have the BC and the RNA count
  BCcount <- data.table(as.data.frame(rev(sort(table(barcodeTable$BC))), row.names = "Var1"), keep.rownames = TRUE)
  setnames(BCcount,colnames(BCcount),c("BC","RNAcount"))
  # set the key to BC
  setkey(BCcount,"BC")
  # extract only BC that are in BCcount
  foundFrags <- output.Table[BCcount,nomatch=0]
  # set the key of foundFrags and LUT.dna to LUTnr
  setkey(foundFrags,"LUTnr")
  setkey(LUT.dna,"LUTnr")
  # extract only the fragments that are in foundFrags and the LUT.dna
  foundFrags <- foundFrags[LUT.dna,nomatch=0]
  # change the name to Sequence and fragement
  setnames(foundFrags,"Sequence","fragment")
  # delete the Name and i.Structure column
  foundFrags[,c("Names","i.Structure"):=NULL]

  # function to to index fragments that are found multiple times since they are the same in multiple sequences from the known Retrograde_transport seq
  matchRange <- function(idxFrag) {
    matchRanges <- which(mcols(allFragments.ranges)$Sequence == foundFrags$fragment[idxFrag])
    return(cbind(matchRanges,idxFrag))
  }

  # create a list of the matches
  match.ranges.list <- mclapply(1:nrow(foundFrags), matchRange, mc.preschedule = TRUE, 
                                mc.cores = detectCores())
  # create a matrix of the matches
  match.ranges <- do.call(rbind, match.ranges.list)
  # create a list of the found fragments
  foundFragments.ranges <- allFragments.ranges[match.ranges[,1]]
  # if there are more than one match, then save the found fragments
  if (ncol(match.ranges) >= 2) {
  foundFrags <- foundFrags[match.ranges[,"idxFrag"],]
  # remove unnecessary columns
  foundFrags[,c("Reads","fragment","Structure","LUTnr"):=NULL]
  # add found fragments to the foundFragments.ranges
  mcols(foundFragments.ranges) <- c(mcols(foundFragments.ranges),foundFrags)
  # sort them by the RNA count in descending order
  foundFragments.ranges <- foundFragments.ranges[order(-mcols(foundFragments.ranges)$RNAcount)]
  # save the found fragments for the sample in the output folder
  saveRDS(foundFragments.ranges, file=paste("2_output/","found.",name.out,".rds", sep=""), 
          compress = TRUE)
  }
  return(log.table)
}

#'Analysis summary
#'============================
# run the analysis for a single tissue
# log.table <- analyzeTissue(2)

# run the analysis for all tissues
log.table <- lapply(1:nrow(load.list), analyzeTissue)

# create a data.table from the log table
log.df <- do.call(rbind, log.table)
# save the log table
saveRDS(log.df, file="S5_log.table.rds", compress = TRUE)

# cleanup the temp files
unlink(paste(tempdir(), "/*", sep = ""), recursive = FALSE, force = FALSE) 

print("Total execution time:")
print(Sys.time()-strt)
devtools::session_info()