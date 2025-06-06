
ggdotplot <- function(hits,
                      type,
                      outDir = NULL,
                      minGenes2plot = 100,
                      appendName = "synHits",
                      dotsPerIn = 256,
                      quantileThresh = .5,
                      plotSize = 12,
                      minScore = 50,
                      maxFacets = 10000,
                      verbose = is.null(outDir)){
  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- blkID <-
    inBuffer <- rnd2 <- rnd1 <- n <- isArrayRep2 <- isArrayRep1 <- chr1 <-
    noAnchor <- bitScore <- quantile <- chr2 <- sameOG <- isAnchor <- NULL
  
  ##############################################################################
  # 1. Get the plot size figured out
  tp <- data.table(hits)
  
  un1 <- uniqueN(tp$ofID1)
  un2 <- uniqueN(tp$ofID2)
  if(un1 > un2){
    ht <- plotSize
    wd <- ht * (un1/un2)
  }else{
    wd <- plotSize
    ht <- wd * (un2/un1)
  }
  
  x <- max(tp$ord1, na.rm = T)
  y <- max(tp$ord2, na.rm = T)
  
  ordPerIn <- x / dotsPerIn
  totDots <- wd * dotsPerIn
  xrnd2 <- floor(x / totDots)+1
  
  ordPerIn <- y / dotsPerIn
  totDots <- ht * dotsPerIn
  yrnd2 <- floor(y / totDots)+1
  
  
  tp[,`:=`(rnd1 = round_toInteger(ord1, xrnd2),
           rnd2 = round_toInteger(ord2, yrnd2))]
  tp <- subset(tp, complete.cases(tp[,c("rnd1", "rnd2", "chr1", "chr2")]))
  
  ##############################################################################
  # 2. Make the plot with all hits, regardless of og
  
  # -- 2.1 subset the hits to those with high enough score
  makeCrappyDotplots <- list(p0 = NULL, p1 = NULL, p2 = NULL)
  if(type %in% c("all", "raw")){
    hc <- subset(tp, bitScore > minScore)
    ng1 <- as.integer(uniqueN(hc$ofID1))
    ng2 <- as.integer(uniqueN(hc$ofID2))
    
    # -- 2.2 get axis labels
    xlab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome1[1], ng1)
    ylab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome2[1], ng2)
    
    # -- 2.3 subset to chrs with enough genes on them
    hc[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1], na.rm = T), by = "chr1"]
    hc[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2], na.rm = T), by = "chr2"]
    hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)
    
    # -- 2.4 count n hits in each aggregated position
    hc <- hc[,c("chr1", "chr2", "rnd1", "rnd2")]
    hc <- subset(hc, complete.cases(hc))
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    setorder(hc, -n)
    hc <- subset(hc, !is.na(n))
    
    # -- 2.5 threshold n to not highlight super strong regions
    qthresh <- quantile(hc$n, quantileThresh)
    if(qthresh > 20)
      qthresh <- 20
    if(qthresh < 5)
      qthresh <- 5
    hc$n[hc$n > qthresh] <- qthresh
    
    # -- 2.6 get plot title
    titlab <- sprintf(
      "All blast hits with score > %s, %s/%s-gene x/y windows (heatmap range: 2-%s+ hits/window)",
      minScore, xrnd2, yrnd2, round(qthresh))
    
    # -- 2.7 make the plot
    setorder(hc, n)
    hc <- subset(hc, n > 1)
    nfacets <- nrow(with(hc, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(tp$chr1[order(tp$rnd1)])
      chrOrd2 <- unique(tp$chr2[order(tp$rnd2)])
      hc[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
               chr2 = factor(chr2, levels = chrOrd2))]
      p0 <- ggplot(hc, aes(rnd1, rnd2, col = n)) +
        geom_point(pch = ".") +
        scale_color_viridis_c(begin = .1, trans = "log10", guide = "none") +
        scale_x_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
        scale_y_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
        theme_genespace()+
        facet_grid(chr2 ~ chr1, scales = "free",
                   space = "free", as.table = F, switch = "both")+
        labs(x = xlab, y = ylab, title = titlab)
    }else{
      p0 <- NULL
      makeCrappyDotplots[["p0"]] <- subset(tp, bitScore > minScore)
    }
    
    ##############################################################################
    # 3. Make the plot with just OG hits
    
    # -- 2.1 subset the hits to those with high enough score
    hc <- subset(tp, sameOG)
    ng1 <- as.integer(uniqueN(hc$ofID1))
    ng2 <- as.integer(uniqueN(hc$ofID2))
    
    # -- 2.2 get axis labels
    xlab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome1[1], ng1)
    ylab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome2[1], ng2)
    
    # -- 2.3 subset to chrs with enough genes on them
    hc[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1]), by = "chr1"]
    hc[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2]), by = "chr2"]
    hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)
    
    # -- 2.4 count n hits in each aggregated position
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    setorder(hc, -n)
    hc <- subset(hc, !is.na(n))
    
    # -- 2.5 threshold n to not highlight super strong regions
    qthresh <- quantile(hc$n, quantileThresh)
    if(qthresh > 20)
      qthresh <- 20
    if(qthresh < 5)
      qthresh <- 5
    hc$n[hc$n > qthresh] <- qthresh
    
    # -- 2.6 get plot title
    titlab <- sprintf(
      "Blast hits where query and target are in the same orthogroup, %s/%s-gene x/y windows (heatmap range: 1-%s+ hits/window)",
      xrnd2, yrnd2, round(qthresh))
    
    # -- 2.7 make the plot
    setorder(hc, n)
    nfacets <- nrow(with(hc, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(tp$chr1[order(tp$rnd1)])
      chrOrd2 <- unique(tp$chr2[order(tp$rnd2)])
      hc[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
               chr2 = factor(chr2, levels = chrOrd2))]
      p1 <- ggplot(hc, aes(rnd1, rnd2, col = n)) +
        geom_point(pch = ".") +
        scale_color_viridis_c(begin = .1, trans = "log10", guide = "none") +
        scale_x_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
        scale_y_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
        theme_genespace()+
        facet_grid(chr2 ~ chr1, scales = "free",
                   space = "free", as.table = F, switch = "both")+
        labs(x = xlab, y = ylab, title = titlab)
    }else{
      p1 <- NULL
      makeCrappyDotplots[["p1"]] <- subset(tp, sameOG)
    }
  }else{
    p1 <- p0 <- NULL
  }
  if(type %in% c("all", "syntenic")){
    ##############################################################################
    # 4. Make the plot with just anchors
    hcBlk <- subset(tp, isAnchor)
    hcBlk[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1]), by = "chr1"]
    hcBlk[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2]), by = "chr2"]
    hcBlk <- subset(hcBlk, ngene1 > minGenes2plot & ngene2 > minGenes2plot)
    hcBlk <- hcBlk[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2", "blkID")]
    blkCols <- sample(gs_colors(uniqueN(hcBlk$blkID)))
    
    ng1 <- as.integer(uniqueN(hcBlk$ofID1))
    ng2 <- as.integer(uniqueN(hcBlk$ofID2))
    
    nfacets <- nrow(with(hcBlk, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(hcBlk$chr1[order(hcBlk$rnd1)])
      chrOrd2 <- unique(hcBlk$chr2[order(hcBlk$rnd2)])
      hcBlk[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
                  chr2 = factor(chr2, levels = chrOrd2))]
      hcBlk <- subset(hcBlk, !is.na(rnd1) & !is.na(rnd2))
      if(nrow(hcBlk) < 1){
        warning(sprintf("no syntenic hits found for %s vs. %s",
                        hits$genome1[1], hits$genome2[1]))
        p2 <- NULL
      }else{
        p2 <- ggplot(hcBlk, aes(x = rnd1, y = rnd2, col = blkID)) +
          geom_point(pch = 15,size=2.2) +
          scale_color_manual(values = blkCols, guide = "none") +
          scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd1), by = 1e3))+
          scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd2), by = 1e3))+
          theme_genespace() +
          facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
          labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                           hits$genome1[1], uniqueN(hits$ofID1[hits$isAnchor])),
               y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                           hits$genome2[1], uniqueN(hits$ofID2[hits$isAnchor])),
               title = sprintf("Syntenic anchor blast hits, colored by block ID"))
      }
    }else{
      p2 <- NULL
      makeCrappyDotplots[["p2"]] <- subset(tp, isAnchor)
    }
  }else{
    p2 <- NULL
  }
  
  if(is.null(outDir)){
    if(verbose)
      cat("writing to the present graphics device")
    if(!is.null(p0))
      print(p0)
    if(!is.null(p1))
      print(p1)
    if(!is.null(p2))
      print(p2)
  }else{
    dpFile <- file.path(outDir,
                        sprintf("%s_vs_%s.%sHits.pdf",
                                tp$genome1[1], tp$genome2[1], type))
    pdf(dpFile, height = ht, width = wd)
    if(verbose)
      cat(sprintf("writing to file: %s", dpFile))
    if(!is.null(makeCrappyDotplots[["p0"]])){
      with(makeCrappyDotplots[["p0"]], plot(
        rnd1, rnd2, pch = ".",
        xlab = "all gene rank order (genome1)",
        ylab = "all gene rank order (genome2)",
        main = "these genomes have too many chromosomes to plot nicely."))
    }else{
      if(!is.null(p0))
        print(p0)
    }
    
    if(!is.null(makeCrappyDotplots[["p0"]])){
      with(makeCrappyDotplots[["p0"]], plot(
        rnd1, rnd2, pch = ".",
        xlab = "orthogroup hit gene rank order (genome1)",
        ylab = "orthogroup gene rank order (genome2)",
        main = "these genomes have too many chromosomes to plot nicely."))
    }else{
      if(!is.null(p0))
        print(p1)
    }
    
    if(!is.null(makeCrappyDotplots[["p2"]])){
      with(makeCrappyDotplots[["p2"]], plot(
        rnd1, rnd2, pch = ".",
        xlab = "syntenic hit gene rank order (genome1)",
        ylab = "syntenic hit gene rank order (genome2)",
        main = "these genomes have too many chromosomes to plot nicely."))
    }else{
      if(!is.null(p2))
        print(p2)
    }
    
    de <- dev.off()
  }
}




theme_genespace <- function(col = "white"){
  theme(panel.background = element_rect(fill = col),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(
          color = "grey50", size = .2, linetype = 2),
        panel.spacing = unit(0, "mm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linewidth = 0.4, colour = "black", fill = NA),
        strip.background = element_blank(),
        strip.text.x = element_text(
          angle = 90, family = "Helvetica", size = 7),
        strip.text.y.left = element_text(
          angle = 0, family = "Helvetica", size = 7),
        axis.title = element_text(family = "Helvetica", size = 8),
        plot.title = element_text(family = "Helvetica", size = 10))
}


gs_colors <- function(n = 10){
  cols <- c("#C4645C", "#F5915C", "#FFC765",
            "#93CD6A", "#7AC7E2", "#66B8FF", "#6666FF", "#9C63E1",
            "#F4BDFF")
  pal <- colorRampPalette(cols)
  return(pal(n))
}
