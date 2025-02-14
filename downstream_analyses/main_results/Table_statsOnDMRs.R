## DMRs statistics from csv files

library(data.table)
library(stringr)
#set correct working directory
setwd("/Volumes/DataStefan/ARPEGGIO_comparisons/LL_syn4Vpro1")
setwd("~/EnvironmentalStressPolyploidEvolution/downstream_analyses/main_results/data/Fig6a/MILD_syn1_v_pro1/dmrseq")

#Read in regions polyploid vs diploid
lyr_kam_CG <- fread("CG_context/parent2_v_allo.txt")
lyr_kam_CHG <- fread("CHG_context/parent2_v_allo.txt")
lyr_kam_CHH <- fread("CHH_context/parent2_v_allo.txt")
print(paste0("The number of significant DMRs in CG context is: ", sum(lyr_kam_CG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHG context is: ", sum(lyr_kam_CHG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHH context is: ", sum(lyr_kam_CHH$qval < 0.05, na.rm = TRUE)))

hal_kam_CG <- fread("CG_context/parent1_v_allo.txt")
hal_kam_CHG <- fread("CHG_context/parent1_v_allo.txt")
hal_kam_CHH <- fread("CHH_context/parent1_v_allo.txt")
print(paste0("The number of significant DMRs in CG context is: ", sum(hal_kam_CG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHG context is: ", sum(hal_kam_CHG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHH context is: ", sum(hal_kam_CHH$qval < 0.05, na.rm = TRUE)))

#Read in regions polyploid vs polyploid
kam_CG <- fread("CG_context/A_v_B_polyploid.txt")
kam_CHG <- fread("CHG_context/A_v_B_polyploid.txt")
kam_CHH <- fread("CHH_context/A_v_B_polyploid.txt")

scaffold_number_CG <- as.integer(str_split_fixed(kam_CG$seqnames, "_", 2)[,2])
scaffold_number_CHG <- as.integer(str_split_fixed(kam_CHG$seqnames, "_", 2)[,2])
scaffold_number_CHH <- as.integer(str_split_fixed(kam_CHH$seqnames, "_", 2)[,2])
hal_kam_CG <- kam_CG[scaffold_number_CG < 2240]
lyr_kam_CG <- kam_CG[scaffold_number_CG > 2239]
hal_kam_CHG <- kam_CHG[scaffold_number_CHG < 2240]
lyr_kam_CHG <- kam_CHG[scaffold_number_CHG > 2240]
hal_kam_CHH <- kam_CHH[scaffold_number_CHH < 2240]
lyr_kam_CHH <- kam_CHH[scaffold_number_CHH > 2240]

print(paste0("The number of significant DMRs in CG context is: ", sum(lyr_kam_CG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHG context is: ", sum(lyr_kam_CHG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHH context is: ", sum(lyr_kam_CHH$qval < 0.05, na.rm = TRUE)))

print(paste0("The number of significant DMRs in CG context is: ", sum(hal_kam_CG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHG context is: ", sum(hal_kam_CHG$qval < 0.05, na.rm = TRUE)))
print(paste0("The number of significant DMRs in CHH context is: ", sum(hal_kam_CHH$qval < 0.05, na.rm = TRUE)))


# Number of significant hyper and hypo DMRs

#lyrata-side
print(paste0("Significant DMRs hypermethylated lyrata-side in CG context: ", 
             sum(lyr_kam_CG$qval < 0.05 & lyr_kam_CG$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypermethylated lyrata-side in CHG context: ", 
             sum(lyr_kam_CHG$qval < 0.05 & lyr_kam_CHG$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypermethylated lyrata-side in CHH context: ", 
             sum(lyr_kam_CHH$qval < 0.05 & lyr_kam_CHH$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated lyrata-side in CG context: ", 
             sum(lyr_kam_CG$qval < 0.05 & lyr_kam_CG$stat > 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated lyrata-side in CHG context: ", 
             sum(lyr_kam_CHG$qval < 0.05 & lyr_kam_CHG$stat > 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated lyrata-side in CHH context: ", 
             sum(lyr_kam_CHH$qval < 0.05 & lyr_kam_CHH$stat > 0, na.rm = TRUE)))

#halleri-side
print(paste0("Significant DMRs hypermethylated halleri-side in CG context: ", 
             sum(hal_kam_CG$qval < 0.05 & hal_kam_CG$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypermethylated halleri-side in CHG context: ", 
             sum(hal_kam_CHG$qval < 0.05 & hal_kam_CHG$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypermethylated halleri-side in CHH context: ", 
             sum(hal_kam_CHH$qval < 0.05 & hal_kam_CHH$stat < 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated halleri-side in CG context: ", 
             sum(hal_kam_CG$qval < 0.05 & hal_kam_CG$stat > 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated halleri-side in CHG context: ", 
             sum(hal_kam_CHG$qval < 0.05 & hal_kam_CHG$stat > 0, na.rm = TRUE)))
print(paste0("Significant DMRs hypomethylated halleri-side in CHH context: ", 
             sum(hal_kam_CHH$qval < 0.05 & hal_kam_CHH$stat > 0, na.rm = TRUE)))

# Select significant hyper and hypo DMRs

# lyrata side
lyr_kam_CG_sig_hyper <- lyr_kam_CG[lyr_kam_CG$qval < 0.05 & lyr_kam_CG$stat < 0, ]
lyr_kam_CG_sig_hypo <- lyr_kam_CG[lyr_kam_CG$qval < 0.05 & lyr_kam_CG$stat > 0, ]
lyr_kam_CG_sig_tot <- lyr_kam_CG[lyr_kam_CG$qval < 0.05, ]

lyr_kam_CHG_sig_hyper <- lyr_kam_CHG[lyr_kam_CHG$qval < 0.05 & lyr_kam_CHG$stat < 0, ]
lyr_kam_CHG_sig_hypo <- lyr_kam_CHG[lyr_kam_CHG$qval < 0.05 & lyr_kam_CHG$stat > 0, ]
lyr_kam_CHG_sig_tot <- lyr_kam_CHG[lyr_kam_CHG$qval < 0.05, ]

lyr_kam_CHH_sig_hyper <- lyr_kam_CHH[lyr_kam_CHH$qval < 0.05 & lyr_kam_CHH$stat < 0, ]
lyr_kam_CHH_sig_hypo <- lyr_kam_CHH[lyr_kam_CHH$qval < 0.05 & lyr_kam_CHH$stat > 0, ]
lyr_kam_CHH_sig_tot <- lyr_kam_CHH[lyr_kam_CHH$qval < 0.05, ]

# halleri side
hal_kam_CG_sig_hyper <- hal_kam_CG[hal_kam_CG$qval < 0.05 & hal_kam_CG$stat < 0, ]
hal_kam_CG_sig_hypo <- hal_kam_CG[hal_kam_CG$qval < 0.05 & hal_kam_CG$stat > 0, ]
hal_kam_CG_sig_tot <- hal_kam_CG[hal_kam_CG$qval < 0.05, ]

hal_kam_CHG_sig_hyper <- hal_kam_CHG[hal_kam_CHG$qval < 0.05 & hal_kam_CHG$stat < 0, ]
hal_kam_CHG_sig_hypo <- hal_kam_CHG[hal_kam_CHG$qval < 0.05 & hal_kam_CHG$stat > 0, ]
hal_kam_CHG_sig_tot <- hal_kam_CHG[hal_kam_CHG$qval < 0.05, ]

hal_kam_CHH_sig_hyper <- hal_kam_CHH[hal_kam_CHH$qval < 0.05 & hal_kam_CHH$stat < 0, ]
hal_kam_CHH_sig_hypo <- hal_kam_CHH[hal_kam_CHH$qval < 0.05 & hal_kam_CHH$stat > 0, ]
hal_kam_CHH_sig_tot <- hal_kam_CHH[hal_kam_CHH$qval < 0.05, ]


#Average length of DMRs

#lyrata
print(paste0("Average length in hyper DMRs CG context (lyrata): ", round(mean(lyr_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CG context (lyrata): ", round(mean(lyr_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CG context (lyrata): ", round(mean(lyr_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Average length in hyper DMRs CHG context (lyrata): ", round(mean(lyr_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CHG context (lyrata): ", round(mean(lyr_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CHG context (lyrata): ", round(mean(lyr_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Average length in hyper DMRs CHH context (lyrata): ", round(mean(lyr_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CHH context (lyrata): ", round(mean(lyr_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CHH context (lyrata): ", round(mean(lyr_kam_CHH_sig_tot$width), 0), "bp"))

#halleri 

print(paste0("Average length in hyper DMRs CG context (halleri): ", round(mean(hal_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CG context (halleri): ", round(mean(hal_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CG context (halleri): ", round(mean(hal_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Average length in hyper DMRs CHG context (halleri): ", round(mean(hal_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CHG context (halleri): ", round(mean(hal_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CHG context (halleri): ", round(mean(hal_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Average length in hyper DMRs CHH context (halleri): ", round(mean(hal_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Average length in hypo DMRs CHH context (halleri): ", round(mean(hal_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Average length in all DMRs CHH context (halleri): ", round(mean(hal_kam_CHH_sig_tot$width), 0), "bp"))


#SD length of DMRs

#lyrata
print(paste0("SD length in hyper DMRs CG context (lyrata): ", round(sd(lyr_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CG context (lyrata): ", round(sd(lyr_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CG context (lyrata): ", round(sd(lyr_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("SD length in hyper DMRs CHG context (lyrata): ", round(sd(lyr_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CHG context (lyrata): ", round(sd(lyr_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CHG context (lyrata): ", round(sd(lyr_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("SD length in hyper DMRs CHH context (lyrata): ", round(sd(lyr_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CHH context (lyrata): ", round(sd(lyr_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CHH context (lyrata): ", round(sd(lyr_kam_CHH_sig_tot$width), 0), "bp"))

#halleri 

print(paste0("SD length in hyper DMRs CG context (halleri): ", round(sd(hal_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CG context (halleri): ", round(sd(hal_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CG context (halleri): ", round(sd(hal_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("SD length in hyper DMRs CHG context (halleri): ", round(sd(hal_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CHG context (halleri): ", round(sd(hal_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CHG context (halleri): ", round(sd(hal_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("SD length in hyper DMRs CHH context (halleri): ", round(sd(hal_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("SD length in hypo DMRs CHH context (halleri): ", round(sd(hal_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("SD length in all DMRs CHH context (halleri): ", round(sd(hal_kam_CHH_sig_tot$width), 0), "bp"))

#Median length of DMRs

#lyrata
print(paste0("Median length in hyper DMRs CG context (lyrata): ", round(median(lyr_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CG context (lyrata): ", round(median(lyr_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CG context (lyrata): ", round(median(lyr_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Median length in hyper DMRs CHG context (lyrata): ", round(median(lyr_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CHG context (lyrata): ", round(median(lyr_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CHG context (lyrata): ", round(median(lyr_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Median length in hyper DMRs CHH context (lyrata): ", round(median(lyr_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CHH context (lyrata): ", round(median(lyr_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CHH context (lyrata): ", round(median(lyr_kam_CHH_sig_tot$width), 0), "bp"))

#halleri 

print(paste0("Median length in hyper DMRs CG context (halleri): ", round(median(hal_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CG context (halleri): ", round(median(hal_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CG context (halleri): ", round(median(hal_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Median length in hyper DMRs CHG context (halleri): ", round(median(hal_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CHG context (halleri): ", round(median(hal_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CHG context (halleri): ", round(median(hal_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Median length in hyper DMRs CHH context (halleri): ", round(median(hal_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Median length in hypo DMRs CHH context (halleri): ", round(median(hal_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Median length in all DMRs CHH context (halleri): ", round(median(hal_kam_CHH_sig_tot$width), 0), "bp"))

#Max length of DMRs

#lyrata
print(paste0("Max length in hyper DMRs CG context (lyrata): ", round(max(lyr_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CG context (lyrata): ", round(max(lyr_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CG context (lyrata): ", round(max(lyr_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Max length in hyper DMRs CHG context (lyrata): ", round(max(lyr_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CHG context (lyrata): ", round(max(lyr_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CHG context (lyrata): ", round(max(lyr_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Max length in hyper DMRs CHH context (lyrata): ", round(max(lyr_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CHH context (lyrata): ", round(max(lyr_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CHH context (lyrata): ", round(max(lyr_kam_CHH_sig_tot$width), 0), "bp"))

#halleri 

print(paste0("Max length in hyper DMRs CG context (halleri): ", round(max(hal_kam_CG_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CG context (halleri): ", round(max(hal_kam_CG_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CG context (halleri): ", round(max(hal_kam_CG_sig_tot$width), 0), "bp"))

print(paste0("Max length in hyper DMRs CHG context (halleri): ", round(max(hal_kam_CHG_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CHG context (halleri): ", round(max(hal_kam_CHG_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CHG context (halleri): ", round(max(hal_kam_CHG_sig_tot$width), 0), "bp"))

print(paste0("Max length in hyper DMRs CHH context (halleri): ", round(max(hal_kam_CHH_sig_hyper$width), 0), "bp"))
print(paste0("Max length in hypo DMRs CHH context (halleri): ", round(max(hal_kam_CHH_sig_hypo$width), 0), "bp"))
print(paste0("Max length in all DMRs CHH context (halleri): ", round(max(hal_kam_CHH_sig_tot$width), 0), "bp"))


#Total length of DMRs

#lyrata
print(paste0("Total length in hyper DMRs CG context (lyrata): ", round(sum(lyr_kam_CG_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CG context (lyrata): ", round(sum(lyr_kam_CG_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CG context (lyrata): ", round(sum(lyr_kam_CG_sig_tot$width)/1000000, 1), "Mb"))

print(paste0("Total length in hyper DMRs CHG context (lyrata): ", round(sum(lyr_kam_CHG_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CHG context (lyrata): ", round(sum(lyr_kam_CHG_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CHG context (lyrata): ", round(sum(lyr_kam_CHG_sig_tot$width)/1000000, 1), "Mb"))

print(paste0("Total length in hyper DMRs CHH context (lyrata): ", round(sum(lyr_kam_CHH_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CHH context (lyrata): ", round(sum(lyr_kam_CHH_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CHH context (lyrata): ", round(sum(lyr_kam_CHH_sig_tot$width)/1000000, 1), "Mb"))

#halleri 

print(paste0("Total length in hyper DMRs CG context (halleri): ", round(sum(hal_kam_CG_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CG context (halleri): ", round(sum(hal_kam_CG_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CG context (halleri): ", round(sum(hal_kam_CG_sig_tot$width)/1000000, 1), "Mb"))

print(paste0("Total length in hyper DMRs CHG context (halleri): ", round(sum(hal_kam_CHG_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CHG context (halleri): ", round(sum(hal_kam_CHG_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CHG context (halleri): ", round(sum(hal_kam_CHG_sig_tot$width)/1000000, 1), "Mb"))

print(paste0("Total length in hyper DMRs CHH context (halleri): ", round(sum(hal_kam_CHH_sig_hyper$width)/1000000, 1), "Mb"))
print(paste0("Total length in hypo DMRs CHH context (halleri): ", round(sum(hal_kam_CHH_sig_hypo$width)/1000000, 1), "Mb"))
print(paste0("Total length in all DMRs CHH context (halleri): ", round(sum(hal_kam_CHH_sig_tot$width)/1000000, 1), "Mb"))


# total Cs included

#lyrata
print(paste0("Total Cs in hyper DMRs CG context (lyrata): ", sum(lyr_kam_CG_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CG context (lyrata): ", sum(lyr_kam_CG_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CG context (lyrata): ", sum(lyr_kam_CG_sig_tot$L)))

print(paste0("Total Cs in hyper DMRs CHG context (lyrata): ", sum(lyr_kam_CHG_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CHG context (lyrata): ", sum(lyr_kam_CHG_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CHG context (lyrata): ", sum(lyr_kam_CHG_sig_tot$L)))

print(paste0("Total Cs in hyper DMRs CHH context (lyrata): ", sum(lyr_kam_CHH_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CHH context (lyrata): ", sum(lyr_kam_CHH_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CHH context (lyrata): ", sum(lyr_kam_CHH_sig_tot$L)))

#halleri
print(paste0("Total Cs in hyper DMRs CG context (halleri): ", sum(hal_kam_CG_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CG context (halleri): ", sum(hal_kam_CG_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CG context (halleri): ", sum(hal_kam_CG_sig_tot$L)))

print(paste0("Total Cs in hyper DMRs CHG context (halleri): ", sum(hal_kam_CHG_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CHG context (halleri): ", sum(hal_kam_CHG_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CHG context (halleri): ", sum(hal_kam_CHG_sig_tot$L)))

print(paste0("Total Cs in hyper DMRs CHH context (halleri): ", sum(hal_kam_CHH_sig_hyper$L)))
print(paste0("Total Cs in hypo DMRs CHH context (halleri): ", sum(hal_kam_CHH_sig_hypo$L)))
print(paste0("Total Cs in all DMRs CHH context (halleri): ", sum(hal_kam_CHH_sig_tot$L)))
