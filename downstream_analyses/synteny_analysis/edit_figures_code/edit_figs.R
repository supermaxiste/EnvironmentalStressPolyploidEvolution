#######################################
### Editing figures for publication ###
#######################################
library(GENESPACE) #GENESPACE v1.3.1: synteny and orthology constrained comparative genomics
library(svglite)
library(ggplot2)
library(data.table)
# Get the snake-GENESPACE working directory:
snake_GENESPACE_dir <- "../snake-GENESPACE"
wd <- paste0(snake_GENESPACE_dir,"/results/genespace/run_dir")
path2mcscanx <- paste0(snake_GENESPACE_dir,"/results/genespace/MCScanX")

# Init genespace: 
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  ploidy = c(1,1,1,1,1),
  genomeIDs = c("A.lyrata_NT1","A.lyrata","A.lyrata_jgi","A.halleri","A.thaliana")
)

### Riparian plots ###
# Change theme for white background:
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
# Plot the figures:
plot_riparian(
  gsParam = gpar,
  useOrder = F,
  refGenome = "A.lyrata_jgi",
  braidAlpha = 0.7,
  chrFill = "lightgrey",
  chrBorderCol = "black",
  addThemes = ggthemes,
  reorderBySynteny = F,
  pdfFile = "../figures/rip_hal_bp.pdf")
plot_riparian(
  gsParam = gpar, 
  refGenome = "A.halleri",
  syntenyWeight = 0,
  braidAlpha = 0.7,
  chrFill = "lightgrey",
  chrBorderCol = "black",
  addThemes = ggthemes,
  pdfFile = "../figures/rip_hal_genOrd.pdf")


### Dotplot ###

# Use modified GENESPACE plotting functions (Make sure it's in your current wd):
source("modified_dotplot.R")
# Read in blast hits:
lyr_v_hal <- read_allBlast(paste0(wd,"/syntenicHits/A.halleri_vs_A.lyrata.allBlast.txt.gz"))
# Make plot:
ggdotplot(lyr_v_hal,type="syntenic",outDir = "../figures/")