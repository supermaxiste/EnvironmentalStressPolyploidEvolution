## This script uses DEG obtained by all the pairwise comparisons
## across diploid and polyploid species in two conditions and plots
## their amount over generations


# Import libraries 

library(tidyverse)
library(patchwork)

# Folder with differentially expressed genes
DEG_folder <- c("path/to/deg_files_folder/")
# Folder with list of low coverage genes
lowC <- c("path/to/downstream_analyses/main_results/data/Fig5a/")

# import low coverage genes list

HM_lowC_genes_hal <- read.table(paste0(lowC, "HM_lowC_genes_hal.txt"), quote="\"", comment.char="")
HM_lowC_genes_lyr <- read.table(paste0(lowC, "HM_lowC_genes_lyr.txt"), quote="\"", comment.char="")

LL_lowC_genes_hal <- read.table(paste0(lowC, "LL_lowC_genes_hal.txt"), quote="\"", comment.char="")
LL_lowC_genes_lyr <- read.table(paste0(lowC, "LL_lowC_genes_lyr.txt"), quote="\"", comment.char="")

# Import files

# HM conditions

HM_halG1_v_G4_DEG <- read.delim(paste0(DEG_folder, "HM_halG1_v_G4_DEG.txt"))
HM_lyrG1_v_G4_DEG <- read.delim(paste0(DEG_folder, "HM_lyrG1_v_G4_DEG.txt"))

HM_halG1_v_RS7G1_DEG <- read.delim(paste0(DEG_folder, "HM_halG1_v_RS7G1_DEG.txt"))
HM_lyrG1_v_RS7G1_DEG <- read.delim(paste0(DEG_folder, "HM_lyrG1_v_RS7G1_DEG.txt"))

HM_halG1_v_RS7G4_DEG <- read.delim(paste0(DEG_folder, "HM_halG1_v_RS7G4_DEG.txt"))
HM_lyrG1_v_RS7G4_DEG <- read.delim(paste0(DEG_folder, "HM_lyrG1_v_RS7G4_DEG.txt"))

HM_RS7G1_v_G4_hal_DEG <- read.delim(paste0(DEG_folder, "HM_RS7G1_v_G4_hal_DEG.txt"))
HM_RS7G1_v_G4_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_RS7G1_v_G4_lyr_DEG.txt"))

HM_ALKG1_v_halG1_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_halG1_DEG.txt"))
HM_ALKG1_v_lyrG1_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_halG1_DEG.txt"))

HM_ALKG1_v_RS7G1_hal_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_RS7G1_hal_DEG.txt"))
HM_ALKG1_v_RS7G1_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_RS7G1_lyr_DEG.txt"))

HM_ALKG1_v_RS7G4_hal_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_RS7G4_hal_DEG.txt"))
HM_ALKG1_v_RS7G4_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_RS7G4_lyr_DEG.txt"))

HM_ALKG1_v_G4_hal_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_G4_hal_DEG.txt"))
HM_ALKG1_v_G4_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_ALKG1_v_G4_lyr_DEG.txt"))

HM_TKSG1_v_halG1_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_halG1_DEG.txt"))
HM_TKSG1_v_lyrG1_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_halG1_DEG.txt"))

HM_TKSG1_v_RS7G1_hal_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_RS7G1_hal_DEG.txt"))
HM_TKSG1_v_RS7G1_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_RS7G1_lyr_DEG.txt"))

HM_TKSG1_v_RS7G4_hal_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_RS7G4_hal_DEG.txt"))
HM_TKSG1_v_RS7G4_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_RS7G4_lyr_DEG.txt"))

HM_TKSG1_v_G5_hal_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_G5_hal_DEG.txt"))
HM_TKSG1_v_G5_lyr_DEG <- read.delim(paste0(DEG_folder, "HM_TKSG1_v_G5_lyr_DEG.txt"))


# LL conditions

LL_halG1_v_G4_DEG <- read.delim(paste0(DEG_folder, "LL_halG1_v_G4_DEG.txt"))
LL_lyrG1_v_G4_DEG <- read.delim(paste0(DEG_folder, "LL_lyrG1_v_G4_DEG.txt"))

LL_halG1_v_RS7G1_DEG <- read.delim(paste0(DEG_folder, "LL_halG1_v_RS7G1_DEG.txt"))
LL_lyrG1_v_RS7G1_DEG <- read.delim(paste0(DEG_folder, "LL_lyrG1_v_RS7G1_DEG.txt"))

LL_halG1_v_RS7G4_DEG <- read.delim(paste0(DEG_folder, "LL_halG1_v_RS7G4_DEG.txt"))
LL_lyrG1_v_RS7G4_DEG <- read.delim(paste0(DEG_folder, "LL_lyrG1_v_RS7G4_DEG.txt"))

LL_RS7G1_v_G4_hal_DEG <- read.delim(paste0(DEG_folder, "LL_RS7G1_v_G4_hal_DEG.txt"))
LL_RS7G1_v_G4_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_RS7G1_v_G4_lyr_DEG.txt"))

LL_ALKG1_v_halG1_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_halG1_DEG.txt"))
LL_ALKG1_v_lyrG1_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_halG1_DEG.txt"))

LL_ALKG1_v_RS7G1_hal_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_RS7G1_hal_DEG.txt"))
LL_ALKG1_v_RS7G1_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_RS7G1_lyr_DEG.txt"))

LL_ALKG1_v_RS7G4_hal_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_RS7G4_hal_DEG.txt"))
LL_ALKG1_v_RS7G4_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_RS7G4_lyr_DEG.txt"))

LL_ALKG1_v_G4_hal_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_G4_hal_DEG.txt"))
LL_ALKG1_v_G4_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_ALKG1_v_G4_lyr_DEG.txt"))

LL_TKSG1_v_halG1_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_halG1_DEG.txt"))
LL_TKSG1_v_lyrG1_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_halG1_DEG.txt"))

LL_TKSG1_v_RS7G1_hal_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_RS7G1_hal_DEG.txt"))
LL_TKSG1_v_RS7G1_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_RS7G1_lyr_DEG.txt"))

LL_TKSG1_v_RS7G4_hal_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_RS7G4_hal_DEG.txt"))
LL_TKSG1_v_RS7G4_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_RS7G4_lyr_DEG.txt"))

LL_TKSG1_v_G5_hal_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_G5_hal_DEG.txt"))
LL_TKSG1_v_G5_lyr_DEG <- read.delim(paste0(DEG_folder, "LL_TKSG1_v_G5_lyr_DEG.txt"))

# We filter low coverage regions since they lead to zero counts in
# fourth generation synthetics and false DEGs

# HM conditions

HM_halG1_v_G4_DEG_flt <- filter(HM_halG1_v_G4_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_lyrG1_v_G4_DEG_flt <- filter(HM_lyrG1_v_G4_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_halG1_v_RS7G1_DEG_flt <- filter(HM_halG1_v_RS7G1_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_lyrG1_v_RS7G1_DEG_flt <- filter(HM_lyrG1_v_RS7G1_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_halG1_v_RS7G4_DEG_flt <- filter(HM_halG1_v_RS7G4_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_lyrG1_v_RS7G4_DEG_flt <- filter(HM_lyrG1_v_RS7G4_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_RS7G1_v_G4_hal_DEG_flt <- filter(HM_RS7G1_v_G4_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_RS7G1_v_G4_lyr_DEG_flt <- filter(HM_RS7G1_v_G4_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_ALKG1_v_halG1_DEG_flt <- filter(HM_ALKG1_v_halG1_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_ALKG1_v_lyrG1_DEG_flt <- filter(HM_ALKG1_v_lyrG1_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_ALKG1_v_RS7G1_hal_DEG_flt <- filter(HM_ALKG1_v_RS7G1_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_ALKG1_v_RS7G1_lyr_DEG_flt <- filter(HM_ALKG1_v_RS7G1_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_ALKG1_v_RS7G4_hal_DEG_flt <- filter(HM_ALKG1_v_RS7G4_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_ALKG1_v_RS7G4_lyr_DEG_flt <- filter(HM_ALKG1_v_RS7G4_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_ALKG1_v_G4_hal_DEG_flt <- filter(HM_ALKG1_v_G4_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_ALKG1_v_G4_lyr_DEG_flt <- filter(HM_ALKG1_v_G4_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_TKSG1_v_halG1_DEG_flt <- filter(HM_TKSG1_v_halG1_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_TKSG1_v_lyrG1_DEG_flt <- filter(HM_TKSG1_v_lyrG1_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_TKSG1_v_RS7G1_hal_DEG_flt <- filter(HM_TKSG1_v_RS7G1_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_TKSG1_v_RS7G1_lyr_DEG_flt <- filter(HM_TKSG1_v_RS7G1_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_TKSG1_v_RS7G4_hal_DEG_flt <- filter(HM_TKSG1_v_RS7G4_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_TKSG1_v_RS7G4_lyr_DEG_flt <- filter(HM_TKSG1_v_RS7G4_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))

HM_TKSG1_v_G5_hal_DEG_flt <- filter(HM_TKSG1_v_G5_hal_DEG, !(geneID %in% HM_lowC_genes_hal$V1))
HM_TKSG1_v_G5_lyr_DEG_flt <- filter(HM_TKSG1_v_G5_lyr_DEG, !(geneID %in% HM_lowC_genes_lyr$V1))


# LL conditions

LL_halG1_v_G4_DEG_flt <- filter(LL_halG1_v_G4_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_lyrG1_v_G4_DEG_flt <- filter(LL_lyrG1_v_G4_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_halG1_v_RS7G1_DEG_flt <- filter(LL_halG1_v_RS7G1_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_lyrG1_v_RS7G1_DEG_flt <- filter(LL_lyrG1_v_RS7G1_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_halG1_v_RS7G4_DEG_flt <- filter(LL_halG1_v_RS7G4_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_lyrG1_v_RS7G4_DEG_flt <- filter(LL_lyrG1_v_RS7G4_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_RS7G1_v_G4_hal_DEG_flt <- filter(LL_RS7G1_v_G4_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_RS7G1_v_G4_lyr_DEG_flt <- filter(LL_RS7G1_v_G4_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_ALKG1_v_halG1_DEG_flt <- filter(LL_ALKG1_v_halG1_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_ALKG1_v_lyrG1_DEG_flt <- filter(LL_ALKG1_v_lyrG1_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_ALKG1_v_RS7G1_hal_DEG_flt <- filter(LL_ALKG1_v_RS7G1_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_ALKG1_v_RS7G1_lyr_DEG_flt <- filter(LL_ALKG1_v_RS7G1_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_ALKG1_v_RS7G4_hal_DEG_flt <- filter(LL_ALKG1_v_RS7G4_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_ALKG1_v_RS7G4_lyr_DEG_flt <- filter(LL_ALKG1_v_RS7G4_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_ALKG1_v_G4_hal_DEG_flt <- filter(LL_ALKG1_v_G4_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_ALKG1_v_G4_lyr_DEG_flt <- filter(LL_ALKG1_v_G4_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_TKSG1_v_halG1_DEG_flt <- filter(LL_TKSG1_v_halG1_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_TKSG1_v_lyrG1_DEG_flt <- filter(LL_TKSG1_v_lyrG1_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_TKSG1_v_RS7G1_hal_DEG_flt <- filter(LL_TKSG1_v_RS7G1_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_TKSG1_v_RS7G1_lyr_DEG_flt <- filter(LL_TKSG1_v_RS7G1_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_TKSG1_v_RS7G4_hal_DEG_flt <- filter(LL_TKSG1_v_RS7G4_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_TKSG1_v_RS7G4_lyr_DEG_flt <- filter(LL_TKSG1_v_RS7G4_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))

LL_TKSG1_v_G5_hal_DEG_flt <- filter(LL_TKSG1_v_G5_hal_DEG, !(geneID %in% LL_lowC_genes_hal$V1))
LL_TKSG1_v_G5_lyr_DEG_flt <- filter(LL_TKSG1_v_G5_lyr_DEG, !(geneID %in% LL_lowC_genes_lyr$V1))


# plotting DEG over generations

# plan for progenitors vs syng1 & g4: stacked barplot with DEG (y-axis) 
# in G1 & G4 (x-axis)

# create dataframe with data from diploids vs synthetics

#HM conditions, halleri side

DEG <- rbind(HM_halG1_v_RS7G1_DEG_flt,
             HM_halG1_v_RS7G4_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("PvsS1", nrow(HM_halG1_v_RS7G1_DEG_flt)),
                                  rep("PvsS4", nrow(HM_halG1_v_RS7G4_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"


HM_gghal <- ggplot(data = DEG, aes(x = generation, fill = as.factor(dir))) +
  geom_bar(width = 0.4) +
  xlab("Generation") +
  ylab("Count") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Cold conditions \nhalleri side")

#HM conditions, lyrata side

DEG_lyr <- rbind(HM_lyrG1_v_RS7G1_DEG_flt,
                 HM_lyrG1_v_RS7G4_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("PvsS1", nrow(HM_lyrG1_v_RS7G1_DEG_flt)),
                                  rep("PvsS4", nrow(HM_lyrG1_v_RS7G4_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

HM_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar(width = 0.4) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")

#LL conditions, halleri side

DEG <- rbind(LL_halG1_v_RS7G1_DEG_flt,
             LL_halG1_v_RS7G4_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("PvsS1", nrow(LL_halG1_v_RS7G1_DEG_flt)),
                                  rep("PvsS4", nrow(LL_halG1_v_RS7G4_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"


LL_gghal <- ggplot(data = DEG, aes(x = generation, fill = as.factor(dir))) +
  geom_bar(width = 0.4) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Hot conditions \nhalleri side")

#LL conditions, lyrata side

DEG_lyr <- rbind(LL_lyrG1_v_RS7G1_DEG_flt,
                 LL_lyrG1_v_RS7G4_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("PvsS1", nrow(LL_lyrG1_v_RS7G1_DEG_flt)),
                                          rep("PvsS4", nrow(LL_lyrG1_v_RS7G4_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

LL_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar(width = 0.4) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")

HM_gghal + HM_gglyr + LL_gghal + LL_gglyr + plot_layout(ncol = 4)

# plan for natural vs all: x-axis shows three time points (pro, syng1, syng4)
# y-axis is number of DEG. Think of ways on how to separate up and down

# create dataframe with data from natural vs all (ALK)

#HM conditions, halleri side

DEG <- rbind(HM_ALKG1_v_halG1_DEG_flt,
             HM_ALKG1_v_RS7G1_hal_DEG_flt,
             HM_ALKG1_v_RS7G4_hal_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("A", nrow(HM_ALKG1_v_halG1_DEG_flt)),
                                  rep("B", nrow(HM_ALKG1_v_RS7G1_hal_DEG_flt)),
                                  rep("C", nrow(HM_ALKG1_v_RS7G4_hal_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"

HM_gghal <- ggplot(data = DEG, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Cold conditions \nhalleri side")


#HM conditions, lyrata side

DEG_lyr <- rbind(HM_ALKG1_v_lyrG1_DEG_flt,
             HM_ALKG1_v_RS7G1_lyr_DEG_flt,
             HM_ALKG1_v_RS7G4_lyr_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("A", nrow(HM_ALKG1_v_lyrG1_DEG_flt)),
                                  rep("B", nrow(HM_ALKG1_v_RS7G1_lyr_DEG_flt)),
                                  rep("C", nrow(HM_ALKG1_v_RS7G4_lyr_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

HM_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")

#LL conditions, halleri side

DEG <- rbind(LL_ALKG1_v_halG1_DEG_flt,
             LL_ALKG1_v_RS7G1_hal_DEG_flt,
             LL_ALKG1_v_RS7G4_hal_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("A", nrow(LL_ALKG1_v_halG1_DEG_flt)),
                                  rep("B", nrow(LL_ALKG1_v_RS7G1_hal_DEG_flt)),
                                  rep("C", nrow(LL_ALKG1_v_RS7G4_hal_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"

LL_gghal <- ggplot(data = DEG, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Hot conditions \nhalleri side")


#LL conditions, lyrata side

DEG_lyr <- rbind(LL_ALKG1_v_lyrG1_DEG_flt,
                 LL_ALKG1_v_RS7G1_lyr_DEG_flt,
                 LL_ALKG1_v_RS7G4_lyr_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("A", nrow(LL_ALKG1_v_lyrG1_DEG_flt)),
                                          rep("B", nrow(LL_ALKG1_v_RS7G1_lyr_DEG_flt)),
                                          rep("C", nrow(LL_ALKG1_v_RS7G4_lyr_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

LL_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("Count") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")


HM_gghal + HM_gglyr + LL_gghal + LL_gglyr + plot_layout(ncol = 4)

# create dataframe with data from natural vs all (TKS)

#HM conditions, halleri side

DEG <- rbind(HM_TKSG1_v_halG1_DEG_flt,
             HM_TKSG1_v_RS7G1_hal_DEG_flt,
             HM_TKSG1_v_RS7G4_hal_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("A", nrow(HM_TKSG1_v_halG1_DEG_flt)),
                                  rep("B", nrow(HM_TKSG1_v_RS7G1_hal_DEG_flt)),
                                  rep("C", nrow(HM_TKSG1_v_RS7G4_hal_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"

HM_gghal <- ggplot(data = DEG, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("Count") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Cold conditions \nhalleri side")


#HM conditions, lyrata side

DEG_lyr <- rbind(HM_TKSG1_v_lyrG1_DEG_flt,
                 HM_TKSG1_v_RS7G1_lyr_DEG_flt,
                 HM_TKSG1_v_RS7G4_lyr_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("A", nrow(HM_TKSG1_v_lyrG1_DEG_flt)),
                                          rep("B", nrow(HM_TKSG1_v_RS7G1_lyr_DEG_flt)),
                                          rep("C", nrow(HM_TKSG1_v_RS7G4_lyr_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

HM_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")

#LL conditions, halleri side

DEG <- rbind(LL_TKSG1_v_halG1_DEG_flt,
             LL_TKSG1_v_RS7G1_hal_DEG_flt,
             LL_TKSG1_v_RS7G4_hal_DEG_flt)

DEG <- mutate(DEG, generation = c(rep("A", nrow(LL_TKSG1_v_halG1_DEG_flt)),
                                  rep("B", nrow(LL_TKSG1_v_RS7G1_hal_DEG_flt)),
                                  rep("C", nrow(LL_TKSG1_v_RS7G4_hal_DEG_flt))))

DEG$dir[DEG$dir == 1] <- "up"
DEG$dir[DEG$dir == -1] <- "down"

LL_gghal <- ggplot(data = DEG, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.position = "none") +
  ylim(0, 12000) +
  ggtitle("Hot conditions \nhalleri side")


#LL conditions, lyrata side

DEG_lyr <- rbind(LL_TKSG1_v_lyrG1_DEG_flt,
                 LL_TKSG1_v_RS7G1_lyr_DEG_flt,
                 LL_TKSG1_v_RS7G4_lyr_DEG_flt)

DEG_lyr <- mutate(DEG_lyr, generation = c(rep("A", nrow(LL_TKSG1_v_lyrG1_DEG_flt)),
                                          rep("B", nrow(LL_TKSG1_v_RS7G1_lyr_DEG_flt)),
                                          rep("C", nrow(LL_TKSG1_v_RS7G4_lyr_DEG_flt))))

DEG_lyr$dir[DEG_lyr$dir == 1] <- "up"
DEG_lyr$dir[DEG_lyr$dir == -1] <- "down"

LL_gglyr <- ggplot(data = DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("NvsP1", "NvsS1", "NvsS4")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ylim(0, 12000) +
  ggtitle("\nlyrata side")


HM_gghal + HM_gglyr + LL_gghal + LL_gglyr + plot_layout(ncol = 4)


# We also plot a barplot with DEGs between G4/5 vs G1

# HM conditions - halleri  side

HM_DEG_hal <- rbind(HM_halG1_v_G4_DEG,
                    HM_RS7G1_v_G4_hal_DEG_flt,
                    HM_ALKG1_v_G4_hal_DEG,
                    HM_TKSG1_v_G5_hal_DEG)

HM_DEG_hal <- mutate(HM_DEG_hal, generation = c(rep("A", nrow(HM_halG1_v_G4_DEG)),
                                                rep("B", nrow(HM_RS7G1_v_G4_hal_DEG_flt)),
                                                rep("C", nrow(HM_ALKG1_v_G4_hal_DEG)),
                                                rep("D", nrow(HM_TKSG1_v_G5_hal_DEG))))

HM_DEG_hal$dir[HM_DEG_hal$dir == 1] <- "up"
HM_DEG_hal$dir[HM_DEG_hal$dir == -1] <- "down"


HM_gghal2 <- ggplot(data = HM_DEG_hal, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("P", "S", "ALK", "TKS")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = "none") +
  ylim(0, 7500) +
  ggtitle("\nhalleri side")

# HM conditions - lyrata side

HM_DEG_lyr <- rbind(HM_lyrG1_v_G4_DEG,
                    HM_RS7G1_v_G4_lyr_DEG_flt,
                    HM_ALKG1_v_G4_lyr_DEG,
                    HM_TKSG1_v_G5_lyr_DEG)

HM_DEG_lyr <- mutate(HM_DEG_lyr, generation = c(rep("A", nrow(HM_lyrG1_v_G4_DEG)),
                                                rep("B", nrow(HM_RS7G1_v_G4_lyr_DEG_flt)),
                                                rep("C", nrow(HM_ALKG1_v_G4_lyr_DEG)),
                                                rep("D", nrow(HM_TKSG1_v_G5_lyr_DEG))))

HM_DEG_lyr$dir[HM_DEG_lyr$dir == 1] <- "up"
HM_DEG_lyr$dir[HM_DEG_lyr$dir == -1] <- "down"


HM_gglyr2 <- ggplot(data = HM_DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("P", "S", "ALK", "TKS")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = "none") +
  ylim(0, 7500) +
  ggtitle("\nlyrata side")

# LL conditions - halleri  side

LL_DEG_hal <- rbind(LL_halG1_v_G4_DEG,
                    LL_RS7G1_v_G4_hal_DEG_flt,
                    LL_ALKG1_v_G4_hal_DEG,
                    LL_TKSG1_v_G5_hal_DEG)

LL_DEG_hal <- mutate(LL_DEG_hal, generation = c(rep("A", nrow(LL_halG1_v_G4_DEG)),
                                                rep("B", nrow(LL_RS7G1_v_G4_hal_DEG_flt)),
                                                rep("C", nrow(LL_ALKG1_v_G4_hal_DEG)),
                                                rep("D", nrow(LL_TKSG1_v_G5_hal_DEG))))

LL_DEG_hal$dir[LL_DEG_hal$dir == 1] <- "up"
LL_DEG_hal$dir[LL_DEG_hal$dir == -1] <- "down"


LL_gghal2 <- ggplot(data = LL_DEG_hal, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("P", "S", "ALK", "TKS")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank(),
        legend.position = "none") +
  ylim(0, 7500) +
  ggtitle("\nhalleri side")

# LL conditions - lyrata side

LL_DEG_lyr <- rbind(LL_lyrG1_v_G4_DEG,
                    LL_RS7G1_v_G4_lyr_DEG_flt,
                    LL_ALKG1_v_G4_lyr_DEG,
                    LL_TKSG1_v_G5_lyr_DEG)

LL_DEG_lyr <- mutate(LL_DEG_lyr, generation = c(rep("A", nrow(LL_lyrG1_v_G4_DEG)),
                                                rep("B", nrow(LL_RS7G1_v_G4_lyr_DEG_flt)),
                                                rep("C", nrow(LL_ALKG1_v_G4_lyr_DEG)),
                                                rep("D", nrow(LL_TKSG1_v_G5_lyr_DEG))))

LL_DEG_lyr$dir[LL_DEG_lyr$dir == 1] <- "up"
LL_DEG_lyr$dir[LL_DEG_lyr$dir == -1] <- "down"


LL_gglyr2 <- ggplot(data = LL_DEG_lyr, aes(x = generation, fill = dir)) +
  geom_bar() +
  scale_x_discrete(labels = c("P", "S", "ALK", "TKS")) +
  xlab("Generation") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size=15),
        legend.title = element_blank()) +
  ylim(0, 7500) +
  ggtitle("\nlyrata side")

(HM_gghal2 + HM_gglyr2) | (LL_gghal2 + LL_gglyr2)
