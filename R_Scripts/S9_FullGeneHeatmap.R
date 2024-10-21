
#' ---
#' title: "Heatmap analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' A secondary approach to visualizing the top candidates and their relative abundance.  
suppressPackageStartupMessages(library(knitr))

#+ setup, include=FALSE, tidy=TRUE
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(kableExtra))
#+ tidy=TRUE
opts_chunk$set(fig.width = 8, fig.height = 10.2)


#'Selection of relevant samples
#'===================
select.samples <- readRDS("data/allSamplesDataTable.RDS")
select.samples$Group[select.samples$Group== "mRNA_3cpc_HEK293T"] <- "mRNA_3cpc_HEK293T"
select.samples$Group[select.samples$Group== "mRNA_30cpc_HEK293T"] <- "mRNA_30cpc_HEK293T"
select.samples <- select.samples[-grep("4wks|mRNA_3cpc_pNeuron_RNA",select.samples$Group),]

select.samples.binCat <- data.table::copy(select.samples)
setkeyv(select.samples.binCat,c("Group","Category"))

select.samples.binCat <- select.samples.binCat[,list(BCcount=length(table(strsplit(paste(t(BC), 
                                                    collapse=","), ","))), NormCount=mean(log2(RNAcount+1))
                                                        ), by=key(select.samples.binCat)]

setkey(select.samples.binCat,Group)
ref.table <- select.samples.binCat["DNA_pscAAVlib"]
ref.table[,c("Group","NormCount"):=NULL]
setnames(ref.table,"BCcount","libBC")
setkey(ref.table,"Category")
select.samples.binCat[,totBC:=sum(BCcount), by="Group"]
max.count <- max(select.samples.binCat$totBC)
select.samples.binCat[,BCcountN:=BCcount/totBC*max.count]
length.Table <- unique(select.samples, by=c("Category","GeneName"))
length.Table <- length.Table[,list(seqlength=sum(seqlength)), by="Category"]
setkey(length.Table,"Category")
setkey(select.samples.binCat,"Category")
select.samples.binCat <- select.samples.binCat[length.Table,nomatch=0]
select.samples.binCat <- select.samples.binCat[ref.table,nomatch=0]

select.samples.binCat[,Category:=gsub("/|_|’","-",Category)]
select.samples.binCat[,BCcountNseq:=BCcountN/seqlength]
select.samples.binCat[,refNormBC:=BCcountN/libBC]


select.samples.binGene <- data.table::copy(select.samples)
setkeyv(select.samples.binGene,c("Group","Category","GeneName","seqlength"))
select.samples.binGene <- select.samples.binGene[,list(BCcount=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                                       NormCount=mean(log2(RNAcount+1))), by=key(select.samples.binGene)]
select.samples.binGene[,GeneName:=gsub("/|_|’","-",GeneName)]
setkey(select.samples.binGene,Group)
ref.table <- select.samples.binGene["DNA_pscAAVlib"]
ref.table[,c("Group","Category","NormCount","seqlength"):=NULL]
setnames(ref.table,"BCcount","libBC")
setkey(ref.table,"GeneName")
setkey(select.samples.binGene,"GeneName")
select.samples.binGene[,totBC:=sum(BCcount), by="Group"]
select.samples.binGene <- select.samples.binGene[ref.table,nomatch=0]
max.count <- max(select.samples.binGene$totBC)
select.samples.binGene[,BCcountN:=BCcount/totBC*max.count]
select.samples.binGene[,BCcountNseq:=BCcountN/seqlength]
select.samples.binGene[,refNormBC:=BCcountN/libBC]

select.samples.binPos <- data.table::copy(select.samples)
setkeyv(select.samples.binPos,c("Group","structure","Sequence"))
select.samples.binPos <- unique(select.samples.binPos, by=c("Group","structure","Sequence")) 
#Due to key, this removes replicates if identical sequence mapped to multiple genes

setkeyv(select.samples.binPos,c("Group","GeneName","AA","seqlength"))
select.samples.binPos <- select.samples.binPos[,list(BCcount=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                                     NormCount=mean(log2(RNAcount+1)),
                                                     AnimalCount=length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                                     LUTnrs=paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                                     mainStruct=paste(unique(structure), collapse=","),
                                                     mismatches=median(mismatches)), by=key(select.samples.binPos)]
setkey(select.samples.binPos,Group)
ref.table <- select.samples.binPos["DNA_pscAAVlib"]
ref.table[,c("Group","NormCount","AnimalCount","LUTnrs","mainStruct","mismatches","seqlength"):=NULL]
setnames(ref.table,"BCcount","libBC")
setkeyv(ref.table,c("GeneName","AA"))
setkeyv(select.samples.binPos,c("GeneName","AA"))
select.samples.binPos <- select.samples.binPos[ref.table,nomatch=0]
select.samples.binPos[,totBC:=sum(BCcount), by="Group"]
max.count <- max(select.samples.binPos$totBC)
select.samples.binPos[,BCcountN:=BCcount/totBC*max.count]
select.samples.binPos[,libNormBC:=BCcountN/libBC]
select.samples.binPos[,BCcountNanim:=BCcountN+AnimalCount]
select.samples.binPos[,BCcountanim:=BCcount+AnimalCount]
select.samples.binPos[,BCcountNseq:=BCcountN/seqlength]
select.samples.binPos[,NormCountBC:=BCcountNseq*NormCount]


select.samples.binPos[,GeneAA:=paste(GeneName," [",AA,"] - ", mainStruct, sep="")]



#'Plot Heatmaps split by Category
#'===================

plotCategory <- function(select.samples.table,plot.col,sample.select){
  setkey(select.samples.table,Group)
  select.samples.select <- select.samples.table[sample.select]
  eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,")", sep="")))
  select.samples.matrix <- acast(select.samples.select, Category~Group, value.var=plot.col) 
  #Utilizes reshape 2 to make matrix for heatmap
  
  select.samples.matrix[is.na(select.samples.matrix)] <- 0
  select.samples.matrix <- select.samples.matrix[,sample.select]
  return(pheatmap(select.samples.matrix, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE))
}

plotCategory(select.samples.binCat,"BCcountNseq",c("DNA_pscAAVlib","mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotCategory(select.samples.binCat,"BCcountNseq",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))
plotCategory(select.samples.binCat,"BCcountNseq",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))

plotCategory(select.samples.binCat,"refNormBC",c("mRNA_30cpc_Str","mRNA_3cpc_Str","mRNA_30cpc_Th","mRNA_3cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_Ctx","mRNA_30cpc_SN","mRNA_3cpc_SN"))
plotCategory(select.samples.binCat,"refNormBC",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))


#'Plot Heatmaps split by GeneName
#'===================

plotGene <- function(select.samples.table,plot.col,sample.select){
  setkey(select.samples.table,Group)
  select.samples.select <- select.samples.table[sample.select]
  eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,")", sep="")))
  select.samples.matrix <- acast(select.samples.select, GeneName~Group, value.var=plot.col) 
  #Utilizes reshape 2 to make matrix for heatmap
  
  select.samples.matrix[is.na(select.samples.matrix)] <- 0
  select.samples.matrix <- select.samples.matrix[,sample.select]
  return(pheatmap(select.samples.matrix, fontsize_row=5.8, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE))
}

plotGene(select.samples.binGene,"BCcountNseq",c("DNA_pscAAVlib","mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotGene(select.samples.binGene,"BCcountNseq",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))
plotGene(select.samples.binGene,"BCcountNseq",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))


plotGene(select.samples.binGene,"refNormBC",c("mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotGene(select.samples.binGene,"refNormBC",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))
plotGene(select.samples.binGene,"refNormBC",c("mRNA_30cpc_Str","mRNA_3cpc_Str","mRNA_30cpc_Th","mRNA_3cpc_Th","mRNA_30cpc_Ctx","mRNA_3cpc_Ctx","mRNA_30cpc_SN","mRNA_3cpc_SN"))
plotGene(select.samples.binGene,"refNormBC",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))

#'Selection of top ten fragments per sample
#'===================
#+ results = 'asis'


setkeyv(select.samples.binPos,c("Group","BCcountNanim","AnimalCount","NormCount"))
setorder(select.samples.binPos,Group,-BCcountanim,-AnimalCount,-BCcount,-NormCount)
setkey(select.samples.binPos,Group)

select.samples.topTwenty <- select.samples.binPos[, head(.SD, 20), by=Group]
select.samples.topTwenty[,c("totBC","seqlength","BCcountN","GeneAA","BCcountNseq"):=NULL]

for (thisGroup in unique(select.samples.topTwenty$Group)){
  out <- select.samples.topTwenty[J(thisGroup)]
  out[,Group:=NULL]
  out[,NormCount:=round(NormCount,digits = 0)]
  setnames(out,c("GeneName","AnimalCount","mismatches","BCcount","NormCount","BCcountNanim"), 
           c(thisGroup,"Animal","missM","BCs","nCount","BCsNan"))
  
  print(knitr::kable(out, format = "latex", booktabs = T) %>%   kable_styling(latex_options = c("striped", "scale_down","repeat_header")) %>% landscape())
}
#+ results = 'markup'


plotPos <- function(select.samples.table,plot.col,sample.select){
  setkeyv(select.samples.table,"Group")
  select.samples.select <- select.samples.table[J(sample.select)]
  eval(parse(text=paste("setorder(select.samples.select,Group, -", plot.col,",-AnimalCount,-NormCount)", sep="")))
  select.samples.topTen <- select.samples.select[, head(.SD, 10), by=Group] 
  select.samples.out <- select.samples.select[select.samples.select$GeneAA %in% select.samples.topTen$GeneAA]
  select.samples.out <- acast(select.samples.out, GeneAA~Group, value.var=plot.col) #Utilizes reshape 2 to make matrix for heatmap
  select.samples.out[is.na(select.samples.out)] <- 0
  select.samples.out <- select.samples.out[,sample.select]
  return(pheatmap(select.samples.out, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE))
}

plotPos(select.samples.binPos,"NormCount",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))
plotPos(select.samples.binPos,"NormCount",c("mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotPos(select.samples.binPos,"NormCount",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))

plotPos(select.samples.binPos,"BCcountNseq",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))
plotPos(select.samples.binPos,"BCcountNseq",c("mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotPos(select.samples.binPos,"BCcountNseq",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))

plotPos(select.samples.binPos,"NormCountBC",c("mRNA_30cpc_pNeuron","mRNA_3cpc_pNeuron","mRNA_30cpc_HEK293T","mRNA_3cpc_HEK293T"))
plotPos(select.samples.binPos,"NormCountBC",c("mRNA_30cpc_Str","mRNA_30cpc_Th","mRNA_30cpc_Ctx","mRNA_30cpc_SN"))
plotPos(select.samples.binPos,"NormCountBC",c("mRNA_3cpc_Str","mRNA_3cpc_Th","mRNA_3cpc_Ctx","mRNA_3cpc_SN"))


devtools::session_info()