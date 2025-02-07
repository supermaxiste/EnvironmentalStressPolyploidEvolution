## Script to plot GO results for our overlaps

# Import libraries

library(tidyverse)
library(ggrepel)

# Import GO results

#GOOverlap_HM_hal_v_RS7G4 <- read.delim("~/OneDrive/PhD/Project/Chapter_3/DMRs_GO_analysis/GOOverlap_HM_hal_v_RS7G4.txt", header=FALSE)
GOOverlap_HM_hal_v_RS7G4 <- read.delim("GOOverlap_HM_hal_v_RS7G4.txt", header=FALSE)
GOOverlap_HM_lyr_v_RS7G1 <- read.delim("GOOverlap_HM_lyr_v_RS7G1.txt", header=FALSE)
GOOverlap_HM_lyr_v_RS7G4 <- read.delim("GOOverlap_HM_lyr_v_RS7G4.txt", header=FALSE)

GOOverlap_LL_hal_v_RS7G1 <- read.delim("GOOverlap_LL_hal_v_RS7G1.txt", header=FALSE)
GOOverlap_LL_hal_v_RS7G4 <- read.delim("GOOverlap_LL_hal_v_RS7G4.txt", header=FALSE)
GOOverlap_LL_lyr_v_RS7G1 <- read.delim("GOOverlap_LL_lyr_v_RS7G1.txt", header=FALSE)
GOOverlap_LL_lyr_v_RS7G4 <- read.delim("GOOverlap_LL_lyr_v_RS7G4.txt", header=FALSE)


# rename columns

rename_columns <- function(overlap_GO){
  
  names(overlap_GO) <- c("GO",
                         "Total",
                         "Overlap",
                         "Expected",
                         "over_under",
                         "FoldEnrichment",
                         "p-value",
                         "FDR")
  return(overlap_GO)
}

GOOverlap_HM_hal_v_RS7G4 <- rename_columns(GOOverlap_HM_hal_v_RS7G4)
GOOverlap_HM_lyr_v_RS7G1 <- rename_columns(GOOverlap_HM_lyr_v_RS7G1)
GOOverlap_HM_lyr_v_RS7G4 <- rename_columns(GOOverlap_HM_lyr_v_RS7G4)

GOOverlap_LL_hal_v_RS7G1 <- rename_columns(GOOverlap_LL_hal_v_RS7G1)
GOOverlap_LL_hal_v_RS7G4 <- rename_columns(GOOverlap_LL_hal_v_RS7G4)
GOOverlap_LL_lyr_v_RS7G1 <- rename_columns(GOOverlap_LL_lyr_v_RS7G1)
GOOverlap_LL_lyr_v_RS7G4 <- rename_columns(GOOverlap_LL_lyr_v_RS7G4)

# GO to label

pointsToLabel_HM_hal_v_RS7G4 <- c("response to inorganic substance (GO:0010035)",
                                  "ion transport (GO:0006811)",
                                  "response to abiotic stimulus (GO:0009628)",
                                  "response to stress (GO:0006950)")

pointsToLabel_HM_lyr_v_RS7G1 <- c("cellular process (GO:0009987)",
                                  "cellular metabolic process (GO:0044237)",
                                  "primary metabolic process (GO:0044238)")

pointsToLabel_HM_lyr_v_RS7G4 <- c("mRNA export from nucleus (GO:0006406)",
                                  "post-embryonic development (GO:0009791)",
                                  "anatomical structure development (GO:0048856)",
                                  "response to stress (GO:0006950)")

pointsToLabel_LL_hal_v_RS7G1 <- c("metabolic process (GO:0008152)",
                                  "organic substance metabolic process (GO:0071704)",
                                  "cellular process (GO:0009987)")

pointsToLabel_LL_hal_v_RS7G4 <- c("RNA splicing (GO:0008380)",
                                  "cell growth (GO:0016049)",
                                  "protein ubiquitination (GO:0016567)",
                                  "cellular process (GO:0009987)")

pointsToLabel_LL_lyr_v_RS7G1 <- c("transport (GO:0006810)",
                                  "metabolic process (GO:0008152)",
                                  "cellular process (GO:0009987)")

pointsToLabel_LL_lyr_v_RS7G4 <- c("metabolic process (GO:0008152)",
                                  "cellular process (GO:0009987)",
                                  "cellular metabolic process (GO:0044237)")

# scatter plot

# HM_hal_v_RS7G4

ggplot(data = GOOverlap_HM_hal_v_RS7G4,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_HM_hal_v_RS7G4[GOOverlap_HM_hal_v_RS7G4$GO %in% pointsToLabel_HM_hal_v_RS7G4,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)

# HM_lyr_v_RS7G1

ggplot(data = GOOverlap_HM_lyr_v_RS7G1,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_HM_lyr_v_RS7G1[GOOverlap_HM_lyr_v_RS7G1$GO %in% pointsToLabel_HM_lyr_v_RS7G1,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)

# HM_lyr_v_RS7G4

ggplot(data = GOOverlap_HM_lyr_v_RS7G4,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_HM_lyr_v_RS7G4[GOOverlap_HM_lyr_v_RS7G4$GO %in% pointsToLabel_HM_lyr_v_RS7G4,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)

# LL_hal_v_RS7G1

ggplot(data = GOOverlap_LL_hal_v_RS7G1,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_LL_hal_v_RS7G1[GOOverlap_LL_hal_v_RS7G1$GO %in% pointsToLabel_LL_hal_v_RS7G1,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)


# LL_hal_v_RS7G4


ggplot(data = GOOverlap_LL_hal_v_RS7G4,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_LL_hal_v_RS7G4[GOOverlap_LL_hal_v_RS7G4$GO %in% pointsToLabel_LL_hal_v_RS7G4,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)

# LL_lyr_v_RS7G1

ggplot(data = GOOverlap_LL_lyr_v_RS7G1,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_LL_lyr_v_RS7G1[GOOverlap_LL_lyr_v_RS7G1$GO %in% pointsToLabel_LL_lyr_v_RS7G1,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)

# LL_lyr_v_RS7G4

ggplot(data = GOOverlap_LL_lyr_v_RS7G4,
       aes(x = FoldEnrichment, y = -log10(FDR))) +
  geom_point(aes(size = Overlap), 
             col = "pink", alpha = 0.5, shape = 16) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_size(range = c(1,20),
             breaks = c(10, 100, 200, 300),
             limits = c(0, 300)) +
  geom_text_repel(data = GOOverlap_LL_lyr_v_RS7G4[GOOverlap_LL_lyr_v_RS7G4$GO %in% pointsToLabel_LL_lyr_v_RS7G4,],
                  aes(label = GO), 
                  box.padding = 1, 
                  size = 6, 
                  max.iter = 1000000, 
                  col = "deeppink") +
  xlim(0,20) +
  ylim(0,15)
