ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTAC
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(beanplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales)) #Gives the log2 ability to ggplot2
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

in.name.P5 <- "data/3em-lib_S3_L004_R1_001.fastq.gz"
in.name.P7 <- "data/3em-lib_S3_L004_R2_001.fastq.gz"
final.name.P5 <- "data/3_S3_L004_R1_001.fastq.gz"
final.name.P7 <- "data/3_S3_L004_R2_001.fastq.gz"

out.name.P5 <- tempfile(pattern = "P5_", tmpdir = tempdir(), fileext = ".fastq.gz")
out.name.P7 <- tempfile(pattern = "P7_", tmpdir = tempdir(), fileext = ".fastq.gz")


command.args <- paste("-Xmx12g overwrite=true k=25 rcomp=f skipr2=t qhdist=0 maskmiddle=t hammingdistance=0 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", shQuote(gsub("([\\])", "", in.name.P5)),
                      " in2=", shQuote(gsub("([\\])", "", in.name.P7)),
                      " outm=", out.name.P5,
                      " outm2=", out.name.P7,
                      " fliteral=", "ATAACTTCGTATAATGTATGCTATACGAATAATTTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGC", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t , ATATCATGGCCGACAAGCAGA

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)

in.name.P5 <- out.name.P5
in.name.P7 <- out.name.P7

colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")

command.args <- paste("-Xmx12g overwrite=true k=10 rcomp=f skipr1=t qhdist=0 maskmiddle=t hammingdistance=0 findbestmatch=t ordered=t threads=",detectCores(),
                      " in=", in.name.P5,
                      " in2=", in.name.P7,
                      " outm=", final.name.P5,
                      " outm2=", final.name.P7,
                      " fliteral=", "CGCCACAACATCGAGGACGGCAGCGTG", sep = "") #Length 48-72 bp k=18 mink=10 qhdist=0 hammingdistance=3 findbestmatch=t , ATATCATGGCCGACAAGCAGA

sys.out <- system2(path.expand("~/bbmap/bbduk2.sh"), args=command.args, stdout=TRUE, stderr=TRUE) #

sys.out <- as.data.frame(sys.out)

colnames(sys.out) <- c("bbduk2 Identification of real amplicons")
invisible(sys.out[" "] <- " ")
lengthOut <- (nrow(sys.out))
knitr::kable(sys.out[3:lengthOut,], format = "markdown")



