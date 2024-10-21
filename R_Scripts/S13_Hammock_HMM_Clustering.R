
#' ---
#' title: "Clustering Polypeptide motifs using the Hammock hidden Markov model peptide clustering"
#' author: "Tomas Bjorklund"
#' header-includes: \usepackage{longtable}
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.6in
#' fontsize: 9.5pt
#' ---

#' This script clusters Polypeptide motifs using the Hammock hidden Markov model peptide clustering.  
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 11, fig.height = 10.5, fig.align = 'center') #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))
strt1<-Sys.time()

#'Loading samples
#'===================

all.samples <- readRDS("data/allSamplesDataTable.RDS")
all.samples[,Peptide:= as.character(Peptide),]

setkey(all.samples,Group)

select.samples <- all.samples[J(c("mRNA_30cpc_SN","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_SN","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_30cpc_SN_4wks","mRNA_30cpc_Th_4wks","mRNA_30cpc_Ctx_4wks","mRNA_3cpc_SN_4wks","mRNA_3cpc_Th_4wks","mRNA_3cpc_Ctx_4wks"))] 

select.samples[,BCcount:=as.integer(mclapply(BC, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Animalcount:=as.integer(mclapply(Animals, function(x) length(table(strsplit(paste(t(x), collapse=","), ","))), mc.cores = detectCores()))]
select.samples[,Score:= BCcount+Animalcount-1,]
select.samples.trsp <- unique(select.samples, by=c("Animals","BC","LUTnrs"))

fasta.names <- paste(1:nrow(select.samples.trsp),select.samples.trsp$Score,select.samples.trsp$Group, sep = "|")
write.fasta(as.list(select.samples.trsp$Peptide), fasta.names, "data/trspSamplesPeptides.fasta", open = "w", nbchar = 60, as.string = TRUE)

#Generate Scoring table for Weblogo Weighting
select.samples.pepMerge <- select.samples.trsp[, sum(Score), by = c("Peptide")]
setnames(select.samples.pepMerge, "V1", "Score")

#'Executing Hammock Clustering
#'===================

Sys.setenv("PATH" = paste("/root/HMMER/binaries",Sys.getenv("PATH"),sep=":"), "HHLIB" = "/home/rstudio/Hammock_v_1.1.1/hhsuite-2.0.16/lib/hh/")
unlink("/home/rstudio/data/Hammock", recursive = TRUE, force = FALSE)
sys.out <-  system(paste("java -jar /home/rstudio/Hammock_v_1.1.1/dist/Hammock.jar full -i /home/rstudio/data/trspSamplesPeptides.fasta -d /home/rstudio/data/Hammock --max_shift 7 -c 65 -t ", detectCores(), sep = ""),
                   intern = TRUE, ignore.stdout = TRUE) 
# Alternative parameters  --use_clinkage --alignment_threshold 23 --max_shift 13 --max_aln_length 37 --count_threshold 50 --max_inner_gaps 0 --assign_thresholds 14.1,10.5,7.0 
hammock.log <- data.table(readLines("data/Hammock/run.log"))

colnames(hammock.log) <- c("Hammock log file")
knitr::kable(hammock.log, longtable = T)


#'Generation of Weblogo visualization
#'===================
ham.clusters <- data.table(read.table("/home/rstudio/data/Hammock/final_clusters.tsv", header = TRUE, skip = 0, sep="\t",
                                      stringsAsFactors = FALSE, fill=TRUE))
id.order <- as.list(ham.clusters$cluster_id)
ham.clusters.all <- data.table(read.table("/home/rstudio/data/Hammock/final_clusters_sequences.tsv", header = TRUE, skip = 0, sep="\t",
                                          stringsAsFactors = FALSE, fill=TRUE))
ham.clusters.all[,alignment := gsub('\\-', '\\_', alignment)]
setkey(select.samples,Peptide)
setkey(select.samples.trsp,Peptide)

unlink("/home/rstudio/data/WEBlogos", recursive = TRUE, force = FALSE)
dir.create(file.path("/home/rstudio/data/", "WEBlogos"), showWarnings = FALSE)
dir.create(file.path("/home/rstudio/data/Hammock/", "alignments_final_Scored"), showWarnings = FALSE)

setkey(ham.clusters.all,cluster_id)
setkey(ham.clusters,cluster_id)
setkey(select.samples.pepMerge,Peptide)


opts_chunk$set(out.width='100%', fig.align = 'center')
generateWeblogo <- function(in.name) {
  #in.name <- ham.clusters$cluster_id[2]
  # in.name <- 6777
  this.fa <- read.fasta(file = paste("/home/rstudio/data/Hammock/alignments_final/", in.name, ".aln", sep=""))
  allSeqs <- unlist(getSequence(this.fa, as.string = TRUE))
  allSeqs <- data.table(unlist(lapply(allSeqs, function(x) gsub("([-])","",toupper(x)))))
  allSeqs.out <- select.samples.pepMerge[J(allSeqs)]
  allSeqs.out$Annot <- data.table(getName(this.fa))
  allSeqs.out[,Annot:= paste(Annot,"_",Score,sep="")]
  allSeqs.out$Alignment <- data.table(toupper(unlist(getSequence(this.fa, as.string = TRUE))))
  
  allSeqs.out <- allSeqs.out[rep(1:.N,Score)][,Indx:=1:.N,by=Peptide]
  allSeqs.out[,Annot:= paste(Annot,"_",Indx,sep="")]
  
  write.fasta(as.list(allSeqs.out$Alignment), allSeqs.out$Annot, nbchar = 60, paste("/home/rstudio/data/Hammock/alignments_final_Scored/", in.name, ".aln", sep=""), open = "w")
  
  
  this.main <- ham.clusters[J(in.name)]
  main.gene <- select.samples.trsp[J(this.main$main_sequence)]$GeneName[1]
  this.title <- paste("## Peptide",this.main$main_sequence,"from",main.gene,"with cluster number",in.name, sep=" ")
  
  tmp <- system(paste("weblogo --format PDF --sequence-type protein --size large --errorbars NO --resolution 299 --composition equiprobable --color-scheme chemistry --title '",this.title,"' < /home/rstudio/data/Hammock/alignments_final_Scored/", in.name, ".aln > /home/rstudio/data/WEBlogos/",in.name,".pdf", sep = ""),
                intern = TRUE, ignore.stdout = FALSE)

  cat('\n')
  cat(this.title, "\n")
  cat('\n')
  cat('\n')
  cat(paste0("![Peptide: ",this.main$main_sequence," from ",main.gene," with cluster number ",in.name,"](/home/rstudio/data/WEBlogos/",in.name,".pdf)"))
  cat('\n')
  out.table <- knitr::kable(this.main[,c(1:8,10)], format = "latex")
  print(column_spec(out.table, 1:9, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  cat('\n')
  this.cluster <- ham.clusters.all[J(in.name)]
  out.table <- knitr::kable(this.cluster[,c(1:7)], format = "latex")
  print(column_spec(out.table, 1:7, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  
  this.found <- select.samples[J(this.cluster$sequence)]
  setnames(this.cluster, "sequence", "Peptide")
  this.found <- merge(this.found,this.cluster[,2:3], by="Peptide", all = FALSE)
  cat('\n')
  cat('\n')
  output.order <- c("alignment","LUTnrs","GeneName","start","structure","Group","Score")
  if (nrow(this.found) >= 48) {
    this.found.p1 <- this.found[1:47,]
    out.table <- knitr::kable(this.found.p1[,..output.order], format = "latex")
    print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "repeat_header")))
    cat('\n')
    cat("\n\n\\pagebreak\n")
    cat("\n\n\\clearpage\n")
    this.found <- this.found[48:nrow(this.found),]
    if (nrow(this.found) >= 48) {
      this.found.p2 <- this.found[1:47,]
      out.table <- knitr::kable(this.found.p2[,..output.order], format = "latex")
      print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "repeat_header")))
      cat('\n')
      cat("\n\n\\pagebreak\n")
      cat("\n\n\\clearpage\n")
      this.found <- this.found[48:nrow(this.found),]
    }
    if (nrow(this.found) >= 48) {
      this.found.p2 <- this.found[1:47,]
      out.table <- knitr::kable(this.found.p2[,..output.order], format = "latex")
      print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "repeat_header")))
      cat('\n')
      cat("\n\n\\pagebreak\n")
      cat("\n\n\\clearpage\n")
      this.found <- this.found[48:nrow(this.found),]
    }
    out.table <- knitr::kable(this.found[,..output.order], format = "latex")
    print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "repeat_header")))
    cat('\n')
    cat("\n\n\\clearpage\n")
    cat('\n')
  } else {
  out.table <- knitr::kable(this.found[,..output.order], format = "latex")
  print(column_spec(out.table, 1, monospace = TRUE) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")))
  cat('\n')
  cat("\n\n\\pagebreak\n")
  cat("\n\n\\clearpage\n")
  }
}

#+ results = 'asis'
invisible(lapply(id.order, generateWeblogo))
#+ results = 'markup'


# setkey(select.samples.trsp,Peptide)
# select.samples.trsp.select <- select.samples.trsp[J(c("PPDELNLTTASLPL"))]



print("Total analysis time:")
print(Sys.time()-strt1)

devtools::session_info()

