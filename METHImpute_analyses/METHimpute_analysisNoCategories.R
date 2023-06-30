## This script will take as input the methylation calls from Bismark and an annotation file
## The output of the script will be the following:
#          1) Conversion rate per context
#          2) Distance correlation plot
#          3) METHimpute model in .RData format
#          4) Number of cytosines: total, covered, M, U, high confidence and low confidence
#          5) Histogram of the methylation proportion for reads with given coverage
#          6) Enrichment plot for covered and missing cytosines

# Import libraries 

library(methimpute)

# Script arguments: first Bismark all C counts (CX_report.txt)
#                   second annotation file in GFF format
#                   third is the prefix for the output files

comm_args <- commandArgs(trailingOnly = TRUE)

# First argument: Bismark all C counts (CX_report.txt)

bismark_allC <- normalizePath(comm_args[1])

# Second argument: annotation file in GFF format

anno <- normalizePath(comm_args[2])

# Third argument: prefix for output files

prefix <- comm_args[3]

# Fourth argument: output folder

output <- normalizePath(comm_args[4])
setwd(output)

# Read Bismark file

allC_data <- importBismark(bismark_allC)
print("Data imported successfully")

# Calculate correlation of methylation levels for all cytosines in all contexts

distcor <- distanceCorrelation(allC_data)
fit <- estimateTransDist(distcor)
print("Correlation computed")

# Save plot as a .png file

corrPlot <- paste0(prefix, "_correlationPlot.png")
png(filename=corrPlot, res = 300, pointsize = 1, width = 8, height = 8, units = "in")
par(cex.axis = 0.1, cex.lab = 0.1)
print(fit$plot)
dev.off()
print("Correlation plots saved")

# Calculate methylation considering context separately 

model <- callMethylationSeparate(data = allC_data, transDist = fit$transDist[c(1,4,6)], num.threads = 8)

# Free some memory

rm(allC_data)

# Get conversion rates and save them

conversion_rates <- 1 - model$params$emissionParams$Unmethylated
write.table(conversion_rates, file = paste0(prefix, "_conversionRates.txt"), quote = FALSE)
print("Conversion rates saved")

# Methylation summary for all cytosines and status

table_all <- table(model$data$status)
write.table(table_all, file = paste0(prefix, "_MU.txt"), quote = FALSE)
table_high_conf <- table(model$data$status[which(model$data$posteriorMax>= 0.9)])
write.table(table_high_conf, file = paste0(prefix, "_high_conf_MU.txt"), quote = FALSE)
print("Methylation summary saved")
status <- cbind(model$data$status, model$data$posteriorMax)
write.table(status, file = paste0(prefix, "_status.txt"), quote = FALSE,
            col.names = FALSE, row.names = FALSE)

# Histogram plotting 

#histPlot <- paste0(prefix, "_histogramPlot.png")
#png(filename=histPlot, res = 300, pointsize = 1, width = 6, height = 6, units = "in")
#par(cex.axis = 0.5, cex.lab = 0.5)
#plotHistogram(model, total.counts = 10)
#dev.off()
#print("Histogram plot saved")

# Enrichment plotting

anno_convert <- toGRanges(anno, format = "GFF")
anno_final <- anno_convert[anno_convert$type=="gene"]

#plot_meth <- plotEnrichment(model$data, annotation=anno_final, range = 2000,
#                            category.column = 'category')

# Instead of plotting we want the underlying dataframe

plot_meth <- plotEnrichment(model$data, annotation=anno_final, range = 500,
                            plot = FALSE)

write.table(plot_meth$meth.lvl, file = paste0(prefix, "_dataframe_meth.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE)

write.table(plot_meth$rc.meth.lvl, file = paste0(prefix, "_dataframe_rcmeth.txt"), quote = FALSE, sep = "\t",
            row.names = FALSE)
# Standard enrichment plot

#enrichmentPlot <- paste0(prefix, "_enrichmentPlot.png")
#png(filename=enrichmentPlot, res = 300, pointsize = 1, width = 10, height = 6, units = "in")
#par(cex.axis = 0.5, cex.lab = 0.5)
#plot_meth$meth.lv
#dev.off()

#print("Enrichment plot 1 saved")

#RCenrichmentPlot <- paste0(prefix, "_RCenrichmentPlot.png")
#png(filename=RCenrichmentPlot, res = 300, pointsize = 1, width = 10, height = 6, units = "in")
#par(cex.axis = 0.5, cex.lab = 0.5)
#plot_meth$rc.meth.lv
#dev.off()

#print("Enrichment plot 2 saved")

# Save image as RData object for potential future use

#save.image(file = paste0(prefix, ".RData"))

#print("RData object saved")
print("The end :)")

