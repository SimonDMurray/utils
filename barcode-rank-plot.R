# This script creates a basic barcode rank plot similar to that seen in cellranger except using Starsolo ouput
# The list of UMIs provided by starsolo in the output file UMIperCellSorted.txt is plotted on the y axis with the row number plotted on the x axis
# Both axes are log scaled
# You can customise plot with colours, annotations and change axis cutoffs etc
library(ggplot2)
umi <- read.table("/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1313/results/FCA_GND8047885/output/GeneFull/UMIperCellSorted.txt")
umi$row_num <- seq.int(nrow(umi))
p <- ggplot(umi, aes(x = umi$row_num, y = umi$V1)) + geom_point() + scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') + xlim(0, 15000) + ylim(0, 15000)
p
