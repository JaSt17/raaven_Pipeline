
#' ---
#' title: "Top 23 connectivity analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This script presents top 23 candidates as connectivity graph.   
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE
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
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(pastecs))
suppressPackageStartupMessages(library(kableExtra))

opts_chunk$set(fig.width = 8, fig.height = 6.5)
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)


#'Selection of relevant samples
#'===================
select.samples <- readRDS("data/allSamplesDataTable.RDS")
select.samples <- select.samples[-grep("4wks|mRNA_3cpc_pNeuron_RNA",select.samples$Group),]

select.samples.merge <- data.table::copy(select.samples)
select.samples.merge$Lib <- "30cpc"
select.samples.merge$Lib[grep("3cpc",select.samples.merge$Group)] <- "3cpc"
select.samples.merge[,Group:=gsub("_30cpc|_3cpc","",Group)]

setkey(select.samples.merge,Group)
select.samples.merge.binPos <- select.samples.merge[c("mRNA_Str","mRNA_Th","mRNA_Ctx","mRNA_SN")]

#'Summarize data over each aa
#'===================
setorder(select.samples.merge.binPos,Group,GeneName,start,width)
winWidth=1
windowTable <- select.samples.merge.binPos[,c("GeneName","start","width"), with = FALSE]
windowTable <- unique(windowTable, by=c("GeneName","start","width"))
windowTable <- windowTable[,seq(start,(start+width-winWidth)),by=c("GeneName","start","width")]
setnames(windowTable,"V1","winStart")
windowTable[,winEnd:=winStart+winWidth-1]
setkeyv(windowTable,c("GeneName","start","width"))
setkeyv(select.samples.merge.binPos,c("GeneName","start","width"))
select.samples.windowBin <- select.samples.merge.binPos[windowTable, nomatch=0, allow.cartesian=TRUE]

setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))
select.samples.windowBin.out <- select.samples.windowBin[, list(Overlaps=.N,
                                                            BCcount=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                                            NormCount=mean(log2(RNAcount+1)),
                                                            AnimalCount=length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                                            LUTnrs=paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                                            mainStruct=paste(unique(structure), collapse=","),
                                                            mainLibs=paste(unique(Lib), collapse=","),
                                                            libCount=length(unique(Lib)),
                                                            mismatches=median(mismatches)
), by=c("Group","GeneName","winStart","winEnd","seqlength")]

#'Identify local maxima for each gene and sample
#'===================

select.samples.windowBin <- data.table::copy(select.samples.windowBin.out)
select.samples.windowBin[,Score:=BCcount+AnimalCount+libCount]
setorder(select.samples.windowBin,Group,-Score,-AnimalCount,-BCcount,-libCount,-NormCount,GeneName,winStart)
select.samples.windowBin[,Rank:=seq(.N), by=c("Group")]
setorder(select.samples.windowBin, Group,GeneName,winStart,-Score)


windowTable <- select.samples.windowBin[,c("Group","GeneName","seqlength"), with = FALSE]
windowTable <- unique(windowTable, by=c("Group","GeneName","seqlength"))
windowTable <- windowTable[,seq(seqlength),by=c("Group","GeneName")]
setnames(windowTable,"V1","winStart")
setkeyv(windowTable,c("Group","GeneName","winStart"))
setkeyv(select.samples.windowBin,c("Group","GeneName","winStart"))
select.samples.windowBin <- select.samples.windowBin[windowTable]
select.samples.windowBin$Score[is.na(select.samples.windowBin$Score)] <- 0

select.samples.windowBin$LUTnrs[is.na(select.samples.windowBin$LUTnrs)] <- seq(length(which(is.na(select.samples.windowBin$LUTnrs))))

setorder(select.samples.windowBin,Group,GeneName,-winStart)
select.samples.windowBin.locMax <- data.table::copy(select.samples.windowBin)
select.samples.windowBin.locMax[,Peak:=as.logical(extract(turnpoints(Score), peak=1, pit=0)),by=c("Group","GeneName")] #,proba=0.0001

setorder(select.samples.windowBin.locMax,Group,GeneName,winStart)
select.samples.windowBin.locMax <- select.samples.windowBin.locMax[select.samples.windowBin.locMax$Peak,]
select.samples.windowBin.locMax[,Peak:=NULL]

setorder(select.samples.windowBin.locMax, Group,Rank,-Score,GeneName,winStart)
select.samples.windowBin.locMax <- unique(select.samples.windowBin.locMax, by=c("Group","LUTnrs"))

#'Make binning for 22aa fragment length
#'===================

#Connect highest scoring fragments

setorder(select.samples.windowBin.locMax,Group,-Score,-AnimalCount,-BCcount,-libCount,-NormCount,GeneName,winStart)
setkey(select.samples.windowBin.locMax,Group)
select.samples.windowBin.locMax[,Rank2:=seq(.N), by=c("Group")]


select.samples.topSelect <- select.samples.windowBin.locMax[, head(.SD, 35), by=Group]


setorder(select.samples.topSelect,Group,GeneName,winStart)
winWidth=21
windowTable <- select.samples.topSelect[,c("GeneName","winStart"), with = FALSE]
windowTable <- unique(windowTable, by=c("GeneName","winStart"))
windowTable <- windowTable[,seq((winStart-winWidth),(winStart)),by=c("GeneName","winStart")]
setnames(windowTable,"V1","binBaseStart")
windowTable <- windowTable[binBaseStart>0,]
windowTable[,binBaseEnd:=binBaseStart+winWidth]
scoreSelect <- select.samples.topSelect[,c("Group","GeneName","winStart","Score"), with = FALSE]
setkeyv(windowTable,c("GeneName","winStart"))
setkeyv(scoreSelect,c("GeneName","winStart"))
scoreSelect <- scoreSelect[windowTable,allow.cartesian=TRUE]

scoreSelect[,mCount:=length(unique(Group)),by=c("GeneName","binBaseStart")]
scoreSelect[,oCount:=.N,by=c("GeneName","binBaseStart")]
scoreSelect[,tCount:=mCount+oCount]
scoreSelect[,offset:=abs(binBaseStart+6-mean(winStart)),c("GeneName","binBaseStart")]
#scoreSelect <- scoreSelect[mCount>oCount,]
setorder(scoreSelect,GeneName,winStart,-tCount,binBaseStart,-offset)

scoreSelect <- scoreSelect[, head(.SD, 1), by=c("GeneName","winStart")]

scoreSelect <- scoreSelect[,c("GeneName","winStart","binBaseStart","binBaseEnd"), with = FALSE]
setkeyv(scoreSelect,c("GeneName","winStart"))
setkeyv(select.samples.windowBin.locMax,c("GeneName","winStart"))


select.samples.windowBin.locMax.bin <- scoreSelect[select.samples.windowBin.locMax,nomatch=0]
setorder(select.samples.windowBin.locMax.bin,GeneName,-binBaseStart,-Score,Group)

select.samples.windowBin.locMerge <- select.samples.windowBin.locMax.bin[, head(.SD, 1), by=c("GeneName","binBaseStart","binBaseEnd","Group")]
select.samples.windowBin.locMerge <- unique(select.samples.windowBin.locMerge, by=c("Group","GeneName","LUTnrs"))

#'Selection of top twenty fragments per sample
#'===================


setorder(select.samples.windowBin.locMerge,Group,-Score,-AnimalCount,-BCcount,-libCount,-NormCount,GeneName,binBaseStart)
setkey(select.samples.windowBin.locMerge,Group)

select.samples.topTwenty <- select.samples.windowBin.locMerge[, head(.SD, 25), by=Group]
select.samples.topTwenty <- select.samples.topTwenty[,c("GeneName","binBaseStart","binBaseEnd"),with=FALSE]
select.samples.topTwenty <- unique(select.samples.topTwenty, by=c("GeneName","binBaseStart"))

setorder(select.samples.windowBin,Group,-Score,-AnimalCount,-BCcount,-libCount,-NormCount,GeneName,winStart)
select.samples.windowBin[,Rank:=seq(.N), by=c("Group")]
setkeyv(scoreSelect,c("GeneName","winStart"))
setkeyv(select.samples.windowBin,c("GeneName","winStart"))

select.samples.windowBin <- scoreSelect[select.samples.windowBin,nomatch=0]
setorder(select.samples.windowBin,Group,GeneName,binBaseStart,-Score)
select.samples.windowBin.Bin <- unique(select.samples.windowBin, by=c("Group","GeneName","binBaseStart"))

setkeyv(select.samples.topTwenty,c("GeneName","binBaseStart"))
setkeyv(select.samples.windowBin,c("GeneName","binBaseStart"))
select.samples.windowBin.allTop <- select.samples.windowBin[select.samples.topTwenty, nomatch=0]
select.samples.windowBin.allTop <- unique(select.samples.windowBin.allTop, by=c("Group","LUTnrs"))
setorder(select.samples.windowBin.allTop,Group,GeneName,binBaseStart,-Score)
select.samples.windowBin.allTop <- unique(select.samples.windowBin.allTop, by=c("Group","GeneName","binBaseStart"))

select.samples.windowBin.allTop[,GeneAA:=paste(GeneName," [",binBaseStart+floor(winWidth/2),"]", sep=""), by=Group]
setorder(select.samples.windowBin.allTop,Group,Rank)

#'Plotting ranked order
#'===================

setkey(select.samples.windowBin.allTop,Group)

v1 <- select.samples.windowBin.allTop["mRNA_Str"]$GeneAA
v2 <- select.samples.windowBin.allTop["mRNA_Th"]$GeneAA
v3 <- select.samples.windowBin.allTop["mRNA_Ctx"]$GeneAA
v4 <- select.samples.windowBin.allTop["mRNA_SN"]$GeneAA


o <- 0.05
DF <- data.table(x = c(rep(1, length(v1)), rep(2, length(v2)),rep(3, length(v3)), rep(4, length(v4))),
                 x1 = c(rep(1, length(v1)),rep(2, length(v2)), rep(3, length(v3)), rep(4, length(v4))),
                 y = c(rev(seq_along(v1)), rev(seq_along(v2)), rev(seq_along(v3)), rev(seq_along(v4))),
                 g = c(v1, v2, v3, v4))
DF[,groupMax:=max(y), by=x]
allMax <- max(DF$y)
DF[,y:=y-groupMax+allMax]


ggplot(DF, aes(x=x, y=y, group=g, label=g)) +
  geom_path(aes(x=x1), 
            size=0.3, color="blue") +
  geom_text(size=1.8) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

