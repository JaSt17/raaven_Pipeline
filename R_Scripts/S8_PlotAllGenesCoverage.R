
#' ---
#' title: "Pairwise sample analysis output"
#' author: "Tomas Bjorklund"
#' edited by: "Jaro Steindorff"
#' output:  
#' ---

#  This will make R-functions such as library() and install.packages() use this directory:
.libPaths(c('~/MyRextensions', .libPaths()))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbio))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(kableExtra))

# start time
strt1<-Sys.time()

#'Generation of infective library
#'===================

# read in the all samples data table
all.samples <- readRDS("0_data/allSamplesDataTable.RDS")

#'Plotting function
#'===================

plotPair <- function(topSample,bottomSample,size.bin=1,winWidth=1,NormalizePlot=TRUE, PlotBC=TRUE) {
  # Select samples
  #===================
  
  # topSample <- "mRNA_3cpc_Ctx"
  # bottomSample <- "mRNA_30cpc_Ctx"
  # filterBC <- FALSE
  # filterAnimal <- FALSE
  # AnimaladjustPlot <- FALSE
  # NormalizePlot <- TRUE
  # size.bin <- 1
  # winWidth=1
  # PlotBC=TRUE
  
  # create the fill values
  fill.values <- eval(parse(text=paste("c(", topSample,"= rgb(38,64,135, maxColorValue = 255), ",
                                       bottomSample,"= rgb(157,190,217, maxColorValue = 255))",sep="")))
  # set the key of all samples to the group
  setkey(all.samples,Group)
  #Select the two samples
  select.samples <- all.samples[Group %in% c(topSample,bottomSample)]
  # adjust the RNA count to log2 and add 1 in case there is a 0
  select.samples[,RNAcount:=log2(RNAcount+1)]
  
  # if the window width is greater than 0
  if (winWidth > 0) {
    # create a ordered table of the samples
    setorder(select.samples,Group,GeneName,start,width)
    #
    windowTable <- select.samples[,c("GeneName","start","width"), with = FALSE]
    windowTable <- unique(windowTable, by=c("GeneName","start","width"))
    windowTable <- windowTable[,(seq(width-winWidth+1)+start-1),by=c("GeneName","start","width")]
    setnames(windowTable,"V1","winStart")
    windowTable[,winEnd:=winStart+winWidth-1]
    setkeyv(windowTable,c("GeneName","start","width"))
    setkeyv(select.samples,c("GeneName","start","width"))
    select.samples.windowBin <- select.samples[windowTable, allow.cartesian=TRUE]
    select.samples.windowBin[,AAproc:=winStart/seqlength*100]
    
    setkey(select.samples.windowBin,Group)
    select.samples.windowBin <- select.samples.windowBin[J(names(fill.values))] #Select the two compared groups
    setkeyv(select.samples.windowBin,c("Group","GeneName","winStart","winEnd"))
    select.samples.windowBin <- select.samples.windowBin[, list(Overlaps=.N,
                                                                seqlength=min(seqlength),
                                                                AAproc = min(AAproc),
                                                                BC = paste(t(BC), collapse=","),
                                                                Animals = paste(t(Animals), collapse=","),
                                                                LUTnrs = paste(t(LUTnrs), collapse=","),
                                                                RNAcount = sum(RNAcount)
    ), by=c("Group","GeneName","winStart","winEnd")]
    
    plot.data.dt <- unique(select.samples.windowBin, by=c("Group","GeneName","winStart","winEnd"))
    
  } 
  else {
    plot.data.dt <- data.table::copy(select.samples)
  }
  

  #===================
  #Binning of data
  #===================
  FullLength <- 100
  position <- seq(0,FullLength,size.bin)
  plot.data.dt[,bin:=findInterval(AAproc, position)]
  
  plot.data.bin <- plot.data.dt[, list(.N,seqlength=min(seqlength),
                                       AAproc = position[findInterval(mean(AAproc),position)],
                                       BCsum=length(table(strsplit(paste(t(BC), collapse=","), ","))),
                                       AnimalCount = length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                                       LUTnrs = paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                                       NormCount = sum(RNAcount)/seqlength*FullLength
  ), by=c("Group","GeneName","bin")]
  plot.data.bin <- unique(plot.data.bin, by=c("Group","GeneName","bin"))
  
  plot.data.bin[,BCanim:=as.double(BCsum*AnimalCount)]
  
  
  #===================
  #Filtration parameters
  #===================
  
  # if the normalization is set to true we will NormCount by dividing by the maximum value
  if (NormalizePlot) {
    plot.data.bin[plot.data.bin$Group == topSample]$NormCount <- plot.data.bin[plot.data.bin$Group == topSample]$NormCount / max(plot.data.bin[plot.data.bin$Group == topSample]$NormCount)
    plot.data.bin[plot.data.bin$Group == bottomSample]$NormCount <- plot.data.bin[plot.data.bin$Group == bottomSample]$NormCount / max(plot.data.bin[plot.data.bin$Group == bottomSample]$NormCount)
  }
  # if the plotBC is set to true we will normalize the BCanim by dividing by the maximum value
  if (PlotBC) {
    plot.data.bin[plot.data.bin$Group == topSample]$BCanim <- plot.data.bin[plot.data.bin$Group == topSample]$BCanim / max(plot.data.bin[plot.data.bin$Group == topSample]$BCanim)
    plot.data.bin[plot.data.bin$Group == bottomSample]$BCanim <- plot.data.bin[plot.data.bin$Group == bottomSample]$BCanim / max(plot.data.bin[plot.data.bin$Group == bottomSample]$BCanim)
  }
  # Flip the values for the bottom samples
  plot.data.bin[Group == bottomSample, `:=`(
    NormCount = NormCount * -1,
    BCanim = BCanim * -1)]


  #===================
  #Output plot
  #===================
  # change the output variable depending on the plotBC
  if (PlotBC) {outVar <- "BCanim"} else {outVar <- "NormCount"}

  # create the ggplot object
  plot.out <- ggplot(plot.data.bin, aes(x = AAproc, y = .data[[outVar]], fill = Group))

  # Add layers
  plot.out <- plot.out + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    scale_fill_manual(name = "Library", values = fill.values) +
    scale_colour_manual(name = "Library", values = fill.values) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    facet_wrap(~ GeneName, ncol = 5) +   
    theme(plot.margin = unit(c(0, 0, 0, 0), units = "mm"),
          legend.position = "bottom",
          legend.spacing = unit(0, "cm"),
          legend.key.height = unit(0, "cm"),
          plot.background = element_rect(fill = "white"),
          axis.text = element_text(size = rel(0.45)),
          axis.ticks = element_line(size = rel(0.5)),
          axis.ticks.length = unit(.05, "cm"),
          strip.text.x = element_text(size = rel(0.5), colour = "black", 
                                      angle = 0, lineheight = 3, vjust = -20),
          strip.background = element_blank(),
          panel.spacing.y = unit(0, "cm"),
          panel.spacing.x = unit(0, "cm"))
  
  #===================
  # Sort and select top samples
  #===================
  # create a copy of the select samples that we will use to sort and select the top samples
  select.samples.binPos <- select.samples
  # set the key to the group, structure and sequence
  setkeyv(select.samples.binPos,c("Group","structure","Sequence"))
  # set the order of the samples
  setorder(select.samples.binPos,Group,structure,Sequence,GeneName)
  # select the unique samples. Due to key, this removes replicates if identical sequence mapped to multiple genes
  select.samples.binPos <- unique(select.samples.binPos, by=c("Group","structure","Sequence")) 
  # set the key to the group, category, gene name and AA
  setkeyv(select.samples.binPos,c("Group","Category","GeneName","AA"))
  # create new columns with the number of BC, the normalized count, the number of animals, the LUT numbers, the main structure and the number of mismatches
  select.samples.binPos[,c("BCcount","NormCount","AnimalCount","LUTnrs","mainStruct","mismatches"):=
                          list(length(table(strsplit(paste(t(BC), collapse=","), ","))),
                               sum(NormCount),
                               length(table(strsplit(paste(t(Animals), collapse=","), ","))),
                               paste(unique(names(table(strsplit(paste(t(LUTnrs), collapse=","), ",")))), collapse=","),
                               paste(unique(structure), collapse=","),
                               median(mismatches)), by=key(select.samples.binPos)]
  # removes duplicate rows, keeping only unique combinations of values in the columns "Group", "NormCount", and "LUTnrs".
  select.samples.binPos <- unique(select.samples.binPos, by=c("Group","NormCount","LUTnrs"))
  # selects and retains only the specified columns and removing any columns not that are not listed.
  select.samples.binPos <- select.samples.binPos[,c("Group","GeneName","AA","NormCount",
                                                    "BCcount","AnimalCount","LUTnrs","mainStruct",
                                                    "mismatches"), with = FALSE]
  
  # creating a order
  # if the plotBC is set to true we will multiply the BC count with the animal count
  if (PlotBC){
    select.samples.binPos[,BCanim:=BCcount*AnimalCount]
    setorder(select.samples.binPos,Group,-BCanim,-BCcount,-AnimalCount,-NormCount)
  } 
  else {
    setorder(select.samples.binPos,Group,-NormCount,-BCcount,-AnimalCount)
  }
  # set the key to group
  setkey(select.samples.binPos,Group)
  # select the top 25 samples
  select.samples.top <- select.samples.binPos[, head(.SD, 25), by=Group]
  topSample <- select.samples.top[J(topSample)]
  bottomSample <- select.samples.top[J(bottomSample)]
  # remove the group column and set the gene name as the column name
  topSample[,c("Group"):=NULL]
  setnames(topSample, "GeneName", topSample)
  bottomSample[,c("Group"):=NULL]
  setnames(bottomSample, "GeneName", bottomSample)
  
  # create the output list with the plot, the bin data and the top and bottom samples
  out.list <- list(plot=plot.out,
                   plotBin=plot.data.bin,
                   top=topSample,
                   bottom=bottomSample)
  
  return(out.list)
}

#'Analyze samples
#'===================

#===================
# Sample plotting
#===================

# function that generates and saves plots for a list of sample pairs
generate_and_save_plots <- function(pairs_list, output_dir = "3_figures", plot_width = 10, plot_height = 10) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop over each pair
  for (pair in pairs_list) {
    # Generate the plot
    out.plot.list <- plotPair(pair[1], pair[2], PlotBC = FALSE)
    plot <- out.plot.list$plot
    
    # Create file name for saving the plot
    file_name <- paste0(output_dir, "/", pair[1], "_vs_", pair[2], "_first.png")
    
    # Save the plot
    ggsave(file_name, plot, width = plot_width, height = plot_height)
  }
}

# comparison of sample pairs
sample_pairs <- list(
  c("DNA_pscAAVlib_Prep2", "DNA_pscAAVlib"),
  c("DNA_AAVlib_DNAse_3cpc", "DNA_AAVlib_DNAse_30cpc"),
  c("DNA_AAVlib_DNAse_30cpc", "DNA_pscAAVlib_Prep2"),
  c("mRNA_30cpc_Trsp", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Th", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Ctx"),
  c("mRNA_3cpc_Trsp", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Th", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Ctx")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# 3cpc vs 30cpc analysis
sample_pairs <- list(
  c("mRNA_30cpc_Trsp", "mRNA_3cpc_Trsp"),
  c("mRNA_30cpc_Str", "mRNA_3cpc_Str"),
  c("mRNA_30cpc_Th", "mRNA_3cpc_Th"),
  c("mRNA_30cpc_Ctx", "mRNA_3cpc_Ctx"),
  c("mRNA_30cpc_SN", "mRNA_3cpc_SN")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# 8wks vs 4wks analysis
sample_pairs <- list(
  c("mRNA_30cpc_Str_4wks", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Th_4wks", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Ctx_4wks", "mRNA_30cpc_Ctx"),
  c("mRNA_30cpc_SN_4wks", "mRNA_30cpc_SN"),
  c("mRNA_3cpc_Str_4wks", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Th_4wks", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Ctx_4wks", "mRNA_3cpc_Ctx"),
  c("mRNA_3cpc_SN_4wks", "mRNA_3cpc_SN")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# In vitro analysis
sample_pairs <- list(
  c("mRNA_3cpc_pNeuron", "mRNA_30cpc_pNeuron"),
  c("mRNA_3cpc_HEK293T", "mRNA_30cpc_HEK293T")
  c("mRNA_All", "DNA_pscAAVlib"),
  c("mRNA_3cpc_Str", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_Th", "mRNA_30cpc_Th"),
  c("mRNA_3cpc_Ctx", "mRNA_30cpc_Ctx"),
  c("mRNA_3cpc_SN", "mRNA_30cpc_SN"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Ctx"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Th", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Ctx"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Th", "mRNA_3cpc_Str"),
  c("mRNA_30cpc_Trsp", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_Trsp", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Trsp", "mRNA_30cpc_Trsp"),
  c("mRNA_3cpc_pNeuron", "mRNA_30cpc_pNeuron"),
  c("mRNA_3cpc_HEK293T", "mRNA_30cpc_HEK293T"),
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# Binning analysis version 2

# function that generates and saves plots for a list of sample pairs but wuth PlotBC set to TRUE
generate_and_save_plots <- function(pairs_list, output_dir = "3_figures", plot_width = 10, plot_height = 10) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop over each pair
  for (pair in pairs_list) {
    # Generate the plot
    out.plot.list <- plotPair(pair[1], pair[2])
    plot <- out.plot.list$plot
    
    # Create file name for saving the plot
    file_name <- paste0(output_dir, "/", pair[1], "_vs_", pair[2], "_second.png")
    
    # Save the plot
    ggsave(file_name, plot, width = plot_width, height = plot_height)
  }
}

# comparison of sample pairs
sample_pairs <- list(
  c("DNA_pscAAVlib_Prep2", "DNA_pscAAVlib"),
  c("DNA_AAVlib_DNAse_3cpc", "DNA_AAVlib_DNAse_30cpc"),
  c("DNA_AAVlib_DNAse_30cpc", "DNA_pscAAVlib_Prep2"),
  c("mRNA_30cpc_Trsp", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Th", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Ctx"),
  c("mRNA_3cpc_Trsp", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Th", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Ctx")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# 3cpc vs 30cpc analysis
sample_pairs <- list(
  c("mRNA_30cpc_Trsp", "mRNA_3cpc_Trsp"),
  c("mRNA_30cpc_Str", "mRNA_3cpc_Str"),
  c("mRNA_30cpc_Th", "mRNA_3cpc_Th"),
  c("mRNA_30cpc_Ctx", "mRNA_3cpc_Ctx"),
  c("mRNA_30cpc_SN", "mRNA_3cpc_SN")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# 8wks vs 4wks analysis
sample_pairs <- list(
  c("mRNA_30cpc_Str_4wks", "mRNA_30cpc_Str"),
  c("mRNA_30cpc_Th_4wks", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Ctx_4wks", "mRNA_30cpc_Ctx"),
  c("mRNA_30cpc_SN_4wks", "mRNA_30cpc_SN"),
  c("mRNA_3cpc_Str_4wks", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Th_4wks", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Ctx_4wks", "mRNA_3cpc_Ctx"),
  c("mRNA_3cpc_SN_4wks", "mRNA_3cpc_SN")
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

# In vitro analysis
sample_pairs <- list(
  c("mRNA_3cpc_pNeuron", "mRNA_30cpc_pNeuron"),
  c("mRNA_3cpc_HEK293T", "mRNA_30cpc_HEK293T")
  c("mRNA_All", "DNA_pscAAVlib"),
  c("mRNA_3cpc_Str", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_Th", "mRNA_30cpc_Th"),
  c("mRNA_3cpc_Ctx", "mRNA_30cpc_Ctx"),
  c("mRNA_3cpc_SN", "mRNA_30cpc_SN"),
  c("mRNA_30cpc_SN", "mRNA_30cpc_Ctx"),
  c("mRNA_30cpc_Ctx", "mRNA_30cpc_Th"),
  c("mRNA_30cpc_Th", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_SN", "mRNA_3cpc_Ctx"),
  c("mRNA_3cpc_Ctx", "mRNA_3cpc_Th"),
  c("mRNA_3cpc_Th", "mRNA_3cpc_Str"),
  c("mRNA_30cpc_Trsp", "mRNA_30cpc_Str"),
  c("mRNA_3cpc_Trsp", "mRNA_3cpc_Str"),
  c("mRNA_3cpc_Trsp", "mRNA_30cpc_Trsp"),
  c("mRNA_3cpc_pNeuron", "mRNA_30cpc_pNeuron"),
  c("mRNA_3cpc_HEK293T", "mRNA_30cpc_HEK293T"),
)
# Generate and save plots
generate_and_save_plots(sample_pairs)

print("Total analysis time:")
print(Sys.time()-strt1)

devtools::session_info()

