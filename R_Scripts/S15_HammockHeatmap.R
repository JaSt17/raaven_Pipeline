
#' ---
#' title: "Heatmaps generated from HMM peptide clustering"
#' author: "Tomas Bjorklund"
#' header-includes: \usepackage{longtable}
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.6in
#' fontsize: 9.5pt
#' ---

#' This script clusters Polypeptide motifs using the Hammock hidden Markov model peptide clustering and generates Heatmaps for most functional motifs.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 11, fig.height = 10.5, fig.align = 'center') #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)



suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
strt1<-Sys.time()

#'Loading samples
#'===================

all.samples <- readRDS("data/allSamplesDataTable.RDS")
all.samples[,Peptide:= as.character(Peptide),]

setkey(all.samples,Group)


#'Generation of heatmaps for in vivo transported samples
#'===================

select.samples <- all.samples[J(c("mRNA_30cpc_Str", "mRNA_30cpc_SN","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_Str", "mRNA_3cpc_SN","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_30cpc_Str_4wks", "mRNA_30cpc_SN_4wks","mRNA_30cpc_Th_4wks","mRNA_30cpc_Ctx_4wks","mRNA_3cpc_Str_4wks","mRNA_3cpc_SN_4wks","mRNA_3cpc_Th_4wks","mRNA_3cpc_Ctx_4wks"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Animalcount:=as.integer(mclapply(Animals, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount+Animalcount-1,]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/invivoSamplesPeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/HammockInVivo", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/invivoSamplesPeptides.fasta -d /home/rstudio/data/HammockInVivo --max_shift 7 -c 250 --alignment_threshold 26 --assign_thresholds 50,40,30 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --max_shift 13 --max_aln_length 37 --count_threshold 50 --max_inner_gaps 0 --assign_thresholds 14.1,10.5,7.0 
hammock.log <- data.table(readLines("data/HammockInVivo/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)


#'Generation of Weblogo visualization
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/HammockInVivo/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))
id.order <- as.list(ham.clusters$cluster_id)
ham.clusters.all <- data.table(read.table("/home/rstudio/data/HammockInVivo/final_clusters_sequences.tsv", header = TRUE, skip = 0, sep="\t",
                                          stringsAsFactors = FALSE, fill=TRUE))
ham.clusters.all[,alignment := gsub('\\-', '\\_', alignment)]
setkey(select.samples,Peptide)
setkey(select.samples.trsp,Peptide)

unlink("/home/rstudio/data/WEBlogosInVivo", recursive = TRUE, force = FALSE)
dir.create(file.path("/home/rstudio/data/", "WEBlogosInVivo"), showWarnings = FALSE)
dir.create(file.path("/home/rstudio/data/HammockInVivo/", "alignments_final_Scored"), showWarnings = FALSE)

setkey(ham.clusters.all,cluster_id)
setkey(ham.clusters,cluster_id)
setkey(select.samples.pepMerge,Peptide)


opts_chunk$set(out.width='100%', fig.align = 'center')
generateWeblogo <- function(in.name) {
  #in.name <- ham.clusters$cluster_id[12]
  # in.name <- 6777
  this.fa <- read.fasta(file = paste("/home/rstudio/data/HammockInVivo/alignments_final/", in.name, ".aln", sep=""))
  allSeqs <- unlist(getSequence(this.fa, as.string = TRUE))
  allSeqs <- data.table(unlist(lapply(allSeqs, function(x) gsub("([-])","",toupper(x)))))
  allSeqs.out <- select.samples.pepMerge[J(allSeqs)]
  allSeqs.out$Annot <- data.table(getName(this.fa))
  allSeqs.out[,Annot:= paste(Annot,"_",Score,sep="")]
  allSeqs.out$Alignment <- data.table(toupper(unlist(getSequence(this.fa, as.string = TRUE))))
  
  allSeqs.out <- allSeqs.out[rep(1:.N,Score)][,Indx:=1:.N,by=Peptide]
  allSeqs.out[,Annot:= paste(Annot,"_",Indx,sep="")]
  
  write.fasta(as.list(allSeqs.out$Alignment), allSeqs.out$Annot, nbchar = 60, paste("/home/rstudio/data/HammockInVivo/alignments_final_Scored/", in.name, ".aln", sep=""), open = "w")
  
  
  this.main <- ham.clusters[J(in.name)]
  main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]
  this.title <- paste("## Peptide",this.main$main_sequence,"from",main.gene,"with cluster number",in.name, sep=" ")
  
  tmp <- system(paste("weblogo --format PDF --sequence-type protein --size large --errorbars NO --resolution 299 --composition equiprobable --color-scheme chemistry --title '",this.title,"' < /home/rstudio/data/HammockInVivo/alignments_final_Scored/", in.name, ".aln > /home/rstudio/data/WEBlogosInVivo/",in.name,".pdf", sep = ""),
                intern = TRUE, ignore.stdout = FALSE)
}

invisible(mclapply(id.order, generateWeblogo, mc.cores = detectCores()))


ham.clusters.merged <- ham.clusters

ham.clusters.merged[,mRNA_Str := mRNA_30cpc_Str + mRNA_3cpc_Str + mRNA_30cpc_Str_4wks + mRNA_3cpc_Str_4wks]
ham.clusters.merged[,mRNA_SN := mRNA_30cpc_SN + mRNA_3cpc_SN + mRNA_30cpc_SN_4wks + mRNA_3cpc_SN_4wks]
ham.clusters.merged[,mRNA_Th := mRNA_30cpc_Th + mRNA_3cpc_Th + mRNA_30cpc_Th_4wks + mRNA_3cpc_Th_4wks]
ham.clusters.merged[,mRNA_Ctx := mRNA_30cpc_Ctx + mRNA_3cpc_Ctx + mRNA_30cpc_Ctx_4wks + mRNA_3cpc_Ctx_4wks]
ham.clusters.merged[,c("mRNA_30cpc_Str", "mRNA_30cpc_SN","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_Str", "mRNA_3cpc_SN","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_30cpc_Str_4wks", "mRNA_30cpc_SN_4wks","mRNA_30cpc_Th_4wks","mRNA_30cpc_Ctx_4wks","mRNA_3cpc_Str_4wks","mRNA_3cpc_SN_4wks","mRNA_3cpc_Th_4wks","mRNA_3cpc_Ctx_4wks") := NULL]


ham.clusters.merged.melt <- melt(ham.clusters.merged, id=c("cluster_id","main_sequence","sum"))
setkeyv(ham.clusters.merged.melt,"variable")
ham.clusters.topTen <- setorder(setDT(ham.clusters.merged.melt), -value)[, head(.SD, 14), keyby = variable]
#ham.clusters.topTen <- ham.clusters.merged.melt[, head(.SD, 15), by=variable] 
ham.clusters.select <- ham.clusters.merged.melt[ham.clusters.merged.melt$cluster_id %in% unique(ham.clusters.topTen$cluster_id)]

ham.clusters.select[, geneName := lapply(main_sequence, function(x) select.samples.trsp[J(x)]$GeneName[1])]
ham.clusters.select[, listName := paste("Pep:",main_sequence,"from",geneName,"cl:",cluster_id, sep=" ")]

select.samples.out <- acast(ham.clusters.select, listName~variable, value.var="value") #Utilizes reshape 2 to make matrix for heatmap
select.samples.out[is.na(select.samples.out)] <- 0
select.samples.out <- select.samples.out[,c(3, 4, 2, 1)]
select.samples.out.scaled <- scale(select.samples.out, center=FALSE, scale=colSums(select.samples.out))
#select.samples.out.scaled <- select.samples.out.scaled[order(round(select.samples.out.scaled[,1],digits = 2),round(select.samples.out.scaled[,2],digits = 2),round(select.samples.out.scaled[,3],digits = 2),round(select.samples.out.scaled[,4],digits = 2),decreasing=TRUE),]
pheatmap(select.samples.out.scaled, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)


#'Generation of heatmaps for in vitro samples
#'===================


select.samples <- all.samples[J(c("mRNA_3cpc_HEK293T", "mRNA_30cpc_HEK293T","mRNA_3cpc_pNeuron","mRNA_30cpc_pNeuron"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Animalcount:=as.integer(mclapply(Animals, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount+Animalcount-1,]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/invitroSamplesPeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/HammockInVitro", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/invitroSamplesPeptides.fasta -d /home/rstudio/data/HammockInVitro --max_shift 7 -c 50 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --max_shift 13 --max_aln_length 37 --count_threshold 50 --max_inner_gaps 0 --assign_thresholds 14.1,10.5,7.0 
hammock.log <- data.table(readLines("data/HammockInVitro/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)


#'Generation of Weblogo visualization
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/HammockInVitro/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))
id.order <- as.list(ham.clusters$cluster_id)
ham.clusters.all <- data.table(read.table("/home/rstudio/data/HammockInVitro/final_clusters_sequences.tsv", header = TRUE, skip = 0, sep="\t",
                                          stringsAsFactors = FALSE, fill=TRUE))
ham.clusters.all[,alignment := gsub('\\-', '\\_', alignment)]
setkey(select.samples,Peptide)
setkey(select.samples.trsp,Peptide)

unlink("/home/rstudio/data/WEBlogosInVitro", recursive = TRUE, force = FALSE)
dir.create(file.path("/home/rstudio/data/", "WEBlogosInVitro"), showWarnings = FALSE)
dir.create(file.path("/home/rstudio/data/HammockInVitro/", "alignments_final_Scored"), showWarnings = FALSE)

setkey(ham.clusters.all,cluster_id)
setkey(ham.clusters,cluster_id)
setkey(select.samples.pepMerge,Peptide)


opts_chunk$set(out.width='100%', fig.align = 'center')
generateWeblogo <- function(in.name) {
  #in.name <- ham.clusters$cluster_id[12]
  # in.name <- 6777
  this.fa <- read.fasta(file = paste("/home/rstudio/data/HammockInVitro/alignments_final/", in.name, ".aln", sep=""))
  allSeqs <- unlist(getSequence(this.fa, as.string = TRUE))
  allSeqs <- data.table(unlist(lapply(allSeqs, function(x) gsub("([-])","",toupper(x)))))
  allSeqs.out <- select.samples.pepMerge[J(allSeqs)]
  allSeqs.out$Annot <- data.table(getName(this.fa))
  allSeqs.out[,Annot:= paste(Annot,"_",Score,sep="")]
  allSeqs.out$Alignment <- data.table(toupper(unlist(getSequence(this.fa, as.string = TRUE))))
  
  allSeqs.out <- allSeqs.out[rep(1:.N,Score)][,Indx:=1:.N,by=Peptide]
  allSeqs.out[,Annot:= paste(Annot,"_",Indx,sep="")]
  
  write.fasta(as.list(allSeqs.out$Alignment), allSeqs.out$Annot, nbchar = 60, paste("/home/rstudio/data/HammockInVitro/alignments_final_Scored/", in.name, ".aln", sep=""), open = "w")
  
  
  this.main <- ham.clusters[J(in.name)]
  main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]
  this.title <- paste("## Peptide",this.main$main_sequence,"from",main.gene,"with cluster number",in.name, sep=" ")
  
  tmp <- system(paste("weblogo --format PDF --sequence-type protein --size large --errorbars NO --resolution 299 --composition equiprobable --color-scheme chemistry --title '",this.title,"' < /home/rstudio/data/HammockInVitro/alignments_final_Scored/", in.name, ".aln > /home/rstudio/data/WEBlogosInVitro/",in.name,".pdf", sep = ""),
                intern = TRUE, ignore.stdout = FALSE)
}

invisible(mclapply(id.order, generateWeblogo, mc.cores = detectCores()))


ham.clusters.merged <- ham.clusters

# ham.clusters.merged[,mRNA_Str := mRNA_30cpc_Str + mRNA_3cpc_Str + mRNA_30cpc_Str_4wks + mRNA_3cpc_Str_4wks]
# ham.clusters.merged[,mRNA_SN := mRNA_30cpc_SN + mRNA_3cpc_SN + mRNA_30cpc_SN_4wks + mRNA_3cpc_SN_4wks]
# ham.clusters.merged[,mRNA_Th := mRNA_30cpc_Th + mRNA_3cpc_Th + mRNA_30cpc_Th_4wks + mRNA_3cpc_Th_4wks]
# ham.clusters.merged[,mRNA_Ctx := mRNA_30cpc_Ctx + mRNA_3cpc_Ctx + mRNA_30cpc_Ctx_4wks + mRNA_3cpc_Ctx_4wks]
# ham.clusters.merged[,c("mRNA_30cpc_Str", "mRNA_30cpc_SN","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_Str", "mRNA_3cpc_SN","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_30cpc_Str_4wks", "mRNA_30cpc_SN_4wks","mRNA_30cpc_Th_4wks","mRNA_30cpc_Ctx_4wks","mRNA_3cpc_Str_4wks","mRNA_3cpc_SN_4wks","mRNA_3cpc_Th_4wks","mRNA_3cpc_Ctx_4wks") := NULL]

library(reshape)
ham.clusters.merged.melt <- melt(ham.clusters.merged, id=c("cluster_id","main_sequence","sum"))
setkeyv(ham.clusters.merged.melt,"variable")
ham.clusters.topTen <- setorder(setDT(ham.clusters.merged.melt), -value)[, head(.SD, 8), keyby = variable]
#ham.clusters.topTen <- ham.clusters.merged.melt[, head(.SD, 15), by=variable] 
ham.clusters.select <- ham.clusters.merged.melt[ham.clusters.merged.melt$cluster_id %in% unique(ham.clusters.topTen$cluster_id)]

ham.clusters.select[, geneName := lapply(main_sequence, function(x) select.samples.trsp[J(x)]$GeneName[1])]
ham.clusters.select[, listName := paste("Pep:",main_sequence,"from",geneName,"cl:",cluster_id, sep=" ")]

select.samples.out <- acast(ham.clusters.select, listName~variable, value.var="value") #Utilizes reshape 2 to make matrix for heatmap
select.samples.out[is.na(select.samples.out)] <- 0
select.samples.out <- select.samples.out[,c(2, 3, 1,4)]
select.samples.out.scaled <- scale(select.samples.out, center=FALSE, scale=colSums(select.samples.out))
#select.samples.out.scaled <- select.samples.out.scaled[order(round(select.samples.out.scaled[,1],digits = 2),round(select.samples.out.scaled[,2],digits = 2),round(select.samples.out.scaled[,3],digits = 2),round(select.samples.out.scaled[,4],digits = 2),decreasing=TRUE),]
pheatmap(select.samples.out.scaled, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)



select.samples <- all.samples[J(c("DNA_pscAAVlib","DNA_pscAAVlib_Prep2","DNA_AAVlib_DNAse_3cpc", "DNA_AAVlib_DNAse_30cpc"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))






#'Clustering DNAse resistant virions
#'===================

select.samples <- all.samples[J(c("DNA_AAVlib_DNAse_3cpc", "DNA_AAVlib_DNAse_30cpc"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/DNAsePeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/HammockDNAse", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/DNAsePeptides.fasta -d /home/rstudio/data/HammockDNAse --max_shift 7 -c 2000 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --alignment_threshold 26 --assign_thresholds 50,40,30
hammock.log <- data.table(readLines("data/HammockDNAse/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)

#'Generation of Scatter plot generation 
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/HammockDNAse/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))



pred.points <- ggplot(data = ham.clusters,
                      aes(x = DNA_AAVlib_DNAse_3cpc,
                          y = DNA_AAVlib_DNAse_30cpc)) +
  labs(x="DNA_AAVlib_DNAse_3cpc",y="DNA_AAVlib_DNAse_30cpc") +
  geom_point()
print(pred.points)



#'Clustering DNAse resistant virions with library
#'===================


select.samples <- all.samples[J(c("DNA_pscAAVlib","DNA_pscAAVlib_Prep2"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/LibDNAsePeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/HammockLibDNAse", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/LibDNAsePeptides.fasta -d /home/rstudio/data/HammockLibDNAse --max_shift 7 -c 2000 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --alignment_threshold 26 --assign_thresholds 50,40,30
hammock.log <- data.table(readLines("data/HammockLibDNAse/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)

#'Generation of Scatter plot generation 
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/HammockLibDNAse/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))



pred.points <- ggplot(data = ham.clusters,
                      aes(x = DNA_pscAAVlib,
                          y = DNA_pscAAVlib_Prep2)) +
  labs(x="DNA_pscAAVlib",y="DNA_pscAAVlib_Prep2") +
  geom_point()
print(pred.points)

#'Clustering DNAse resistant virions with library
#'===================


select.samples <- all.samples[J(c("DNA_pscAAVlib_Prep2","DNA_AAVlib_DNAse_3cpc","DNA_AAVlib_DNAse_30cpc"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/LibDNAsePeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/HammockLibDNAse", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/LibDNAsePeptides.fasta -d /home/rstudio/data/HammockLibDNAse --max_shift 7 -c 2000 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --alignment_threshold 26 --assign_thresholds 50,40,30
hammock.log <- data.table(readLines("data/HammockLibDNAse/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)

#'Generation of Scatter plot generation 
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/HammockLibDNAse/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))



pred.points <- ggplot(data = ham.clusters,
                      aes(x = DNA_pscAAVlib_Prep2,
                          y = DNA_AAVlib_DNAse_30cpc)) +
  labs(x="DNA_pscAAVlib_Prep2",y="DNA_AAVlib_DNAse_30cpc") +
  geom_point()
print(pred.points)


print("Total analysis time:")
print(Sys.time()-strt1)

devtools::session_info()