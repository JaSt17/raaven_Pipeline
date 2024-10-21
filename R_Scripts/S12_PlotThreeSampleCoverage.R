
#' ---
#' title: "Three sample pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' output: 
#'  pdf_document:
#'    highlight: tango
#' geometry: margin=0.7in
#' fontsize: 10pt
#' ---

#' This is the final script presenting top candidates and overview plots. 
suppressPackageStartupMessages(library(knitr))
#+ setup, include=FALSE

opts_chunk$set(fig.width = 7.3, fig.height = 10.1) #Full height 11
opts_chunk$set(tidy=TRUE)
opts_chunk$set(comment = NA)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggplus))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

#'Generation of infective library
#'===================
all.samples <- readRDS("data/allSamplesDataTable.RDS")

all.samples$Group[all.samples$Group== "mRNA_3cpc_HEK293T"] <- "mRNA_3cpc_HEK293T"
all.samples$Group[all.samples$Group== "mRNA_30cpc_HEK293T"] <- "mRNA_30cpc_HEK293T"

#'Plotting function
#'===================

# Select samples
#===================

  topSampleOne <- "mRNA_3cpc_Th"
  topSampleTwo <- "mRNA_3cpc_Ctx"
  topSampleThree <- "mRNA_3cpc_SN"
  bottomSampleOne <- "mRNA_30cpc_Th"
  bottomSampleTwo <- "mRNA_30cpc_Ctx"
  bottomSampleThree <- "mRNA_30cpc_SN"
  filterBC <- FALSE
  filterAnimal <- FALSE
  AnimaladjustPlot <- FALSE
  NormalizePlot <- TRUE
  size.bin <- 1
  winWidth=1
  PlotBC=TRUE
  
fill.values <- eval(parse(text=paste("c(", topSampleOne,"= rgb(93,52,27, maxColorValue = 255), ",
                                     topSampleTwo,"= rgb(155,98,60, maxColorValue = 255), ",
                                     topSampleThree,"= rgb(213,168,132, maxColorValue = 255), ",
                                     bottomSampleOne,"= rgb(38,64,135, maxColorValue = 255),",
                                     bottomSampleTwo,"= rgb(75,132,165, maxColorValue = 255),",
                                     bottomSampleThree,"= rgb(164,189,216, maxColorValue = 255))",sep="")))
setkey(all.samples,Group)
select.samples <- all.samples[J(names(fill.values))] #Select the six compared groups
select.samples[,RNAcount:=log2(RNAcount+1)]

  setorder(select.samples,Group,GeneName,start,width)
  
  windowTable <- select.samples[,c("GeneName","start","width"), with = FALSE]
  windowTable <- unique(windowTable, by=c("GeneName","start","width"))
  windowTable <- windowTable[,(seq(start,start+width-winWidth)),by=c("GeneName","start","width")]
  setnames(windowTable,"V1","winStart")
  windowTable[,winEnd:=winStart+winWidth-1]
  setkeyv(windowTable,c("GeneName","start","width"))
  setkeyv(select.samples,c("GeneName","start","width"))
  select.samples.windowBin <- select.samples[windowTable, allow.cartesian=TRUE]

setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))
select.samples.windowBin <- select.samples.windowBin[, list(Overlaps=.N,
                                        BC = paste(t(BC), collapse=","),
                                        Animals = paste(t(Animals), collapse=","),
                                        LUTnrs = paste(t(LUTnrs), collapse=","),
                                        RNAcount = sum(RNAcount)
), by=c("Group","GeneName","winStart","winEnd","seqlength")]

plot.data.dt <- unique(select.samples.windowBin, by=c("Group","GeneName","winStart","winEnd"))

#===================
#Binning of data
#===================
FullLength <- max(plot.data.dt$winStart)
position <- seq(0,FullLength,size.bin)
plot.data.dt[,bin:=findInterval(winStart, position)]

plot.data.bin <- plot.data.dt[, list(.N,seqlength=min(seqlength),
                                     BCsum=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                     AA = position[findInterval(mean(winStart),position)],
                                     AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                     LUTnrs = paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                     NormCount = sum(RNAcount)
), by=c("Group","GeneName","bin")]
plot.data.bin <- unique(plot.data.bin, by=c("Group","GeneName","bin"))

plot.data.bin[,BCanim:=as.double(BCsum+AnimalCount)]



#===================
#Filtration parameters
#===================

if (NormalizePlot) {
  for (this.name in names(fill.values)){
  plot.data.bin[plot.data.bin$Group == this.name]$NormCount <- plot.data.bin[plot.data.bin$Group == this.name]$NormCount / max(plot.data.bin[plot.data.bin$Group == this.name]$NormCount)
  }
  }


if (PlotBC && NormalizePlot) {
  for (this.name in names(fill.values)){
  plot.data.bin[plot.data.bin$Group == this.name]$BCanim <- plot.data.bin[plot.data.bin$Group == this.name]$BCanim / max(plot.data.bin[plot.data.bin$Group == this.name]$BCanim)
  }
}

for (this.name in names(fill.values)[seq((length(fill.values)/2)+1,length(fill.values))]){
plot.data.bin[plot.data.bin$Group == this.name]$NormCount <- plot.data.bin[plot.data.bin$Group == this.name]$NormCount*-1 #This line flips the values for the second half of the groups
plot.data.bin[plot.data.bin$Group == this.name]$BCanim <- plot.data.bin[plot.data.bin$Group == this.name]$BCanim*-1 
}
#===================
#Output plot
#===================

if (PlotBC) {outVar <- "BCanim"} else {outVar <- "NormCount"}

plot.out <- eval(parse(text=paste("ggplot(plot.data.bin,aes(x=AA,y=",outVar,", fill = Group))",sep="")))
plot.out <- plot.out + 
  geom_bar(stat="identity", position="identity")+theme_bw() +
  scale_fill_manual(name = "Library", values = fill.values) +
  scale_colour_manual(name = "Library", values = fill.values) +
  scale_x_continuous(breaks=c(seq(0,3000,100)),expand =c(0,0)) +   
  theme(plot.margin=unit(x=c(0,0,0,0),units="mm"),
        legend.position="bottom",
        legend.spacing=unit(0,"cm"),
        legend.key.height=unit(0, "cm"),
        plot.background=element_rect(fill="white"),
        axis.text = element_text(size = rel(0.45)),
        axis.ticks = element_line(size = rel(0.5)),
        axis.ticks.length = unit(.05, "cm"),
        strip.text.x = element_text(size = rel(0.5), colour = "black", 
                                    angle = 0, lineheight=0.1, vjust=0.1),
        strip.background = element_blank(),
        panel.spacing.y = unit(-0.15, "cm"),
        panel.spacing.x = unit(0, "cm"))

facet_multiple(plot = plot.out, facets = 'GeneName', ncol = 1, nrow = 25)

devtools::session_info()

