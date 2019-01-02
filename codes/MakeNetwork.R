###
#   File name : MakeNetwork.R
#   Author    : Hyunjin Kim
#   Date      : Dec 30, 2017
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate Spearman corrleation and build a network based on the correlation
#
#   Instruction
#               1. Source("MakeNetwork.R")
#               2. Run the function "makeCorNet" - specify the input files (microbe & antibiotic)
#               3. The network results will be generated in RedeR container
#
#   Example
#               > source("The_directory_of_MakeNetwork.R/MakeNetwork.R")
#               > makeCorNet(dataPath="./data/argtestre_re.csv", threshold=0.8)
###

makeCorNet <- function(dataPath="./data/argtestre_re.csv", threshold=0.8) {
  
  ### load library
  if(!require(RedeR)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("RedeR")
    library(RedeR)
  }
  
  
  ### load datasets
  d <- read.csv(dataPath, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### keep type of the rows
  t <- d$Type
  names(t) <- rownames(d)
  d <- d[,-1]
  
  
  ### build a correlation matrix with Spearman
  corMat <- cor(t(d), method = "spearman", use = "complete.obs")
  corMat[lower.tri(corMat, diag = TRUE)] <- 0
  corMat[which(corMat < threshold, arr.ind = TRUE)] <- 0
  
  
  ### compute average expressions for the rows
  expAvg <- apply(d, 1, mean)
  
  
  ### A function to remove zero degree nodes
  removeZeroDegreeNodes<-function(mat) {
    zeroDegree <- 0
    cnt <- 1
    for(i in 1:nrow(mat)) {
      if((sum(mat[i,], na.rm = TRUE) == 0) && (sum(mat[,i], na.rm = TRUE) == 0)) {
        zeroDegree[cnt] <- i
        cnt <- cnt + 1
      }
    }
    
    return (mat[-zeroDegree, -zeroDegree])
  }
  
  ### remove zero degree nodes
  corMat <- removeZeroDegreeNodes(corMat)
  t2 <- t[rownames(corMat)]
  expAvg <- expAvg[rownames(corMat)]
  
  
  ### set igraph features
  g <- graph.adjacency(corMat, mode = "upper", weighted = TRUE)
  E(g)$width <- E(g)$weight*5
  E(g)$edgeColor <- "gray"
  g$legEdgeColor$scale <- "gray"
  #V(g)$nodeSize <- (expAvg / max(expAvg)) * 30
  #V(g)$legNodeSize <- V(g)$nodeSize
  V(g)$nodeLineColor <- "black"
  V(g)[which(t2 == "microbe")]$color <- "red"
  V(g)[which(t2 == "antibiotic")]$color <- "green"
  V(g)$legNodeColor <- V(g)$color
  V(g)$nodeFontSize <- as.integer(rep(20, nrow(corMat)))

  
  ### create network
  rdp<-RedPort()
  calld(rdp)
  
  ### add features to the palette
  addGraph(rdp,g, layout.kamada.kawai(g))
  colors<-colorRampPalette(colors=c("red","green"))(2)
  addLegend.color(rdp, colors, labvec=c("Microbe", "Antibiotic"))
  #circleLabel<-floor(seq(min(expAvg),max(expAvg),(max(expAvg) - min(expAvg))/4))
  #circleSize<-(circleLabel / max(circleLabel)) * 30
  #addLegend.size(rdp,sizevec=circleSize,labvec=circleLabel,title="Average Expression")
  shape <- c("ELLIPSE", "DIAMOND")
  addLegend.shape(rdp, shape, title=sprintf("Correlation Network (Spearman > %s)", threshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
  relax(rdp, p3=350, p4=300, p8=5, ps=T)
  
  
  ### HRSD, CA, DEMON & EBBR
  d_hrsd <- d[,1:4]
  d_ca <- d[,9:10]
  d_de <- d[,5:8]
  
  
  ### HRSD
  #
  ### build a correlation matrix with Spearman
  corMat <- cor(t(d_hrsd), method = "spearman", use = "complete.obs")
  corMat[lower.tri(corMat, diag = TRUE)] <- 0
  corMat[which(corMat < threshold, arr.ind = TRUE)] <- 0
  #
  ### remove zero degree nodes
  corMat <- removeZeroDegreeNodes(corMat)
  t2 <- t[rownames(corMat)]
  #
  ### set igraph features
  g <- graph.adjacency(corMat, mode = "upper", weighted = TRUE)
  E(g)$width <- E(g)$weight*5
  E(g)$edgeColor <- "gray"
  g$legEdgeColor$scale <- "gray"
  V(g)$nodeLineColor <- "black"
  V(g)[which(t2 == "microbe")]$color <- "red"
  V(g)[which(t2 == "antibiotic")]$color <- "green"
  V(g)$legNodeColor <- V(g)$color
  V(g)$nodeFontSize <- as.integer(rep(20, nrow(corMat)))
  #
  ### create network
  rdp<-RedPort()
  calld(rdp)
  #
  ### add features to the palette
  addGraph(rdp,g, layout.kamada.kawai(g))
  addLegend.color(rdp, colors, labvec=c("Microbe", "Antibiotic"))
  addLegend.shape(rdp, shape, title=sprintf("HRSD Correlation Network (Spearman > %s)", threshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
  relax(rdp, p3=350, p4=300, p8=5, ps=T)
  
  
  ### CA
  #
  ### build a correlation matrix with Spearman
  corMat <- cor(t(d_ca), method = "spearman", use = "complete.obs")
  corMat[lower.tri(corMat, diag = TRUE)] <- 0
  corMat[which(corMat < threshold, arr.ind = TRUE)] <- 0
  #
  ### remove zero degree nodes
  corMat <- removeZeroDegreeNodes(corMat)
  t2 <- t[rownames(corMat)]
  #
  ### set igraph features
  g <- graph.adjacency(corMat, mode = "upper", weighted = TRUE)
  E(g)$width <- E(g)$weight*5
  E(g)$edgeColor <- "gray"
  g$legEdgeColor$scale <- "gray"
  V(g)$nodeLineColor <- "black"
  V(g)[which(t2 == "microbe")]$color <- "red"
  V(g)[which(t2 == "antibiotic")]$color <- "green"
  V(g)$legNodeColor <- V(g)$color
  V(g)$nodeFontSize <- as.integer(rep(20, nrow(corMat)))
  #
  ### create network
  rdp<-RedPort()
  calld(rdp)
  #
  ### add features to the palette
  addGraph(rdp,g, layout.kamada.kawai(g))
  addLegend.color(rdp, colors, labvec=c("Microbe", "Antibiotic"))
  addLegend.shape(rdp, shape, title=sprintf("CA Correlation Network (Spearman > %s)", threshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
  relax(rdp, p3=350, p4=300, p8=5, ps=T)
  
  
  ### DEMON & EBBR
  #
  ### build a correlation matrix with Spearman
  corMat <- cor(t(d_de), method = "spearman", use = "complete.obs")
  corMat[lower.tri(corMat, diag = TRUE)] <- 0
  corMat[which(corMat < threshold, arr.ind = TRUE)] <- 0
  #
  ### remove zero degree nodes
  corMat <- removeZeroDegreeNodes(corMat)
  t2 <- t[rownames(corMat)]
  #
  ### set igraph features
  g <- graph.adjacency(corMat, mode = "upper", weighted = TRUE)
  E(g)$width <- E(g)$weight*5
  E(g)$edgeColor <- "gray"
  g$legEdgeColor$scale <- "gray"
  V(g)$nodeLineColor <- "black"
  V(g)[which(t2 == "microbe")]$color <- "red"
  V(g)[which(t2 == "antibiotic")]$color <- "green"
  V(g)$legNodeColor <- V(g)$color
  V(g)$nodeFontSize <- as.integer(rep(20, nrow(corMat)))
  #
  ### create network
  rdp<-RedPort()
  calld(rdp)
  #
  ### add features to the palette
  addGraph(rdp,g, layout.kamada.kawai(g))
  addLegend.color(rdp, colors, labvec=c("Microbe", "Antibiotic"))
  addLegend.shape(rdp, shape, title=sprintf("DEMON & EBBR Correlation Network (Spearman > %s)", threshold), position="topleft", ftsize=25, vertical=FALSE, dxtitle=50, dxborder=10, dyborder=-60)
  relax(rdp, p3=350, p4=300, p8=5, ps=T)
  
}


