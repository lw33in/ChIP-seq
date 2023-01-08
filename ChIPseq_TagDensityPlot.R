# Create ChIP-seq tag density plot using ggplot 2

#=====================================================================================
# Directories 
#=====================================================================================
datadir = "/opt/projects/../ChIPseq/"
resultdir = "/opt/projects/../ChIPseq/analysis/"

#=====================================================================================
# Load packages
#=====================================================================================
library(datasets)
library(ggplot2)

#=====================================================================================
# Load tag files 
#=====================================================================================
setwd(datadir)
tag = read.delim("TF_ChIPseq_tag_normalized.txt", header=T)

#=====================================================================================
# Create tag density plots
#=====================================================================================
setwd(resultdir)
plot(tag)
pdf("./Figure1_TF_TagDensityPlot.pdf",width=4.5,height=4,paper='special')
op = par(mar = c(3, 3, 2, 2)) # set graphic margins
plot.new()
plot.window(xlim = c(-2000, 2000), ylim = c(0, 18)) # plot window - x and y axes scales 
# add axes
axis(side = 1, pos = 0, at = seq(from = -2000, to = 2000, by = 1000), col = "black", 
     lwd.ticks = 0.3, cex.axis = 1, col.axis = "black", lwd = 1)
axis(side = 2, pos = "left", at = seq(from = 0, to = 18, by = 6), col = "black", 
     las = 2, lwd.ticks = 0.3, cex.axis =1, col.axis = "black", lwd = 1)
polygon(data_density_plot$distance, data_density_plot$tag, col = "#FF8000",  lwd = 4, border = "#FF1F50" )
# add legends
text(1.3, 0.35, labels = "Normal(0, 1)", col = "gray30", cex = 0.9)
text(3.8, 0.15, labels = "Normal(1, 1.5)", col = "gray30", cex = 0.9)
text(-3.6, 0.15, labels = "Normal(-1, 1.5)", col = "gray30", cex = 0.9)
# add title 
title(main="TF",
      xlab="DistancefromTSS", 
      ylab="ChIP-FragmentDepth")
dev.off()
par(op) # turn off par

