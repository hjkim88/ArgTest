###
#   File name : CorNetCirPlot.R
#   Author    : Hyunjin Kim
#   Date      : Jan 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Make a circular plot based on correlations using circos
#
#   Instruction
#               1. Source("CorNetCirPlot.R")
#               2. Run the function "cirPlot" - specify the input files and output path
#               3. The circular plot will be generated in the output path
#
#   Example
#               > source("The_directory_of_CorNetCirPlot.R/CorNetCirPlot.R")
#               > cirPlot(exp, nodeType, corMethod="spearman", threshold=0.8, outputPath="./results/correlation_network_circular_plot.pdf")
###


cirPlot <- function(exp, nodeType, corMethod="spearman", threshold=0.8, outputPath="./results/correlation_network_circular_plot.pdf") {
  
  ### load library
  if(!require(OmicCircos)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("OmicCircos")
    library(OmicCircos)
  }
  
  
  ### build a correlation matrix
  corMat <- cor(t(exp), method = corMethod, use = "complete.obs")
  
  ### set useless cell values to zero (use only upper triangle)
  corMat[lower.tri(corMat, diag = TRUE)] <- 0
  corMat[which(corMat < threshold, arr.ind = TRUE)] <- 0
  
  
  # ### A function to remove zero degree nodes
  # removeZeroDegreeNodes<-function(mat) {
  #   zeroDegree <- NULL
  #   cnt <- 1
  #   for(i in 1:nrow(mat)) {
  #     if((sum(mat[i,], na.rm = TRUE) == 0) && (sum(mat[,i], na.rm = TRUE) == 0)) {
  #       zeroDegree[cnt] <- i
  #       cnt <- cnt + 1
  #     }
  #   }
  #   
  #   if(length(zeroDegree) > 0) {
  #     return(mat[-zeroDegree, -zeroDegree])
  #   } else {
  #     return(mat)
  #   }
  # }
  # 
  # 
  # ### remove zero degree nodes
  # corMat <- removeZeroDegreeNodes(corMat)
  # nodeType <- nodeType[rownames(corMat)]
  
  
  ### set parameters for circos
  seg.num <- nrow(corMat)
  sample.num <- ncol(exp)
  seg.name <- rownames(corMat)
  fileName <- strsplit(outputPath, "/", fixed = TRUE)
  fileName <- fileName[[1]][length(fileName[[1]])]
  fileName <- strsplit(fileName, ".", fixed = TRUE)[[1]][1]
  
  ### set seg.f
  seg.f <- matrix("NA",seg.num*sample.num, 5)
  colnames(seg.f) <- c("seg.name", "seg.start", "seg.end", "the.v", "NO")
  
  for(i in 1:seg.num) {
    for(j in 1:sample.num) {
      seg.f[(i-1)*sample.num+j,1] <- seg.name[i]
      seg.f[(i-1)*sample.num+j,2] <- j-1
      seg.f[(i-1)*sample.num+j,3] <- j
    }
  }
  
  seg.f <- data.frame(seg.f)
  
  ### set seg.v
  seg.v <- matrix(0,seg.num*sample.num, 3)
  colnames(seg.v) <- c("seg.name", "sample", "exp")
  
  for(i in 1:seg.num) {
    for(j in 1:sample.num) {
      seg.v[(i-1)*sample.num+j,1] <- seg.name[i]
      seg.v[(i-1)*sample.num+j,2] <- j
      seg.v[(i-1)*sample.num+j,3] <- exp[rownames(corMat)[i],j]
    }
  }
  
  seg.v <- data.frame(seg.v)
  
  
  ### set db
  db <- segAnglePo(seg.f, seg=seg.name)
  
  ### set link.pg.v
  link.pg.v <- data.frame(matrix(0, length(which(corMat >= threshold)), 6))
  colnames(link.pg.v) <- c("seg1", "start1", "end1", "seg2", "start2", "end2")
  
  cnt <- 1
  for(i in 1:seg.num) {
    for(j in i:seg.num) {
      if(corMat[i,j] >= threshold) {
        link.pg.v[cnt,1] <- seg.name[i]
        link.pg.v[cnt,2] <- 1+(sample.num-1)/3
        link.pg.v[cnt,3] <- 1+2*(sample.num-1)/3
        link.pg.v[cnt,4] <- seg.name[j]
        link.pg.v[cnt,5] <- 1+(sample.num-1)/3
        link.pg.v[cnt,6] <- 1+2*(sample.num-1)/3
        cnt <- cnt+1
      }
    }
  }
  
  
  ### draw circular plot and save as pdf
  pdf(outputPath, 9, 9)
  par(mar=c(2, 2, 2, 2))
  line_colors <- rainbow(nrow(link.pg.v), alpha=0.3)
  node_colors <- colorRampPalette(colors=c("red","green"))(length(unique(nodeType)))
  names(node_colors) <- unique(nodeType)
  colorType <- node_colors[nodeType]
  plot(c(1,900), c(1,900), type="n", axes=FALSE, xlab="", ylab="", main=fileName)
  circos(xc=450, yc=390, R=300, cir=db, type="chr", col=colorType, print.chr.lab=TRUE, W=4, scale=FALSE)
  circos(xc=450, yc=390, R=260, cir=db, W=40, mapping=seg.v, col.v=3, type="l", B=TRUE, col="black", lwd=2, scale=FALSE)
  circos(xc=450, yc=390, R=250, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(line_colors,nrow(link.pg.v), replace = TRUE))
  dev.off()
  
}



