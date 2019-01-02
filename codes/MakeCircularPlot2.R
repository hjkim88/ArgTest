###
#   File name : MakeCircularPlot2.R
#   Author    : Hyunjin Kim
#   Date      : Apr 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Draw a circular plot based on the data
#
#   Instruction
#               1. Source("MakeCircularPlot2.R")
#               2. Run the function "makeCirPlot2" - specify the input files and output directory
#               3. The circular plot will be generated in the output directory
#
#   Example
#               > source("The_directory_of_MakeCircularPlot2.R/MakeCircularPlot2.R")
#               > makeCirPlot2(dataDir="./data/", outputDir="./results/")
###

makeCirPlot2 <- function(dataDir="./data/", outputDir="./results/") {
  
  ### load library
  if(!require(OmicCircos)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("OmicCircos")
    library(OmicCircos)
  }
  
  
  ### load datasets
  d1 <- read.csv(paste0(dataDir, "data1.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  d2 <- read.csv(paste0(dataDir, "data2.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  d3 <- read.csv(paste0(dataDir, "data3.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### make symmetric matrix
  makeSymmetricMat <- function(df) {
    num <- nrow(df) + ncol(df)
    
    temp1 <- data.frame(matrix(0, ncol(df), ncol(df)))
    colnames(temp1) <- colnames(df)
    rownames(temp1) <- colnames(df)
    
    temp2 <- data.frame(matrix(0, num, nrow(df)))
    colnames(temp2) <- rownames(df)
    rownames(temp2)[1:nrow(df)] <- rownames(df)
    rownames(temp2)[(nrow(df)+1):num] <- colnames(df)
    
    new_df <- rbind(df, temp1)
    new_df <- cbind(temp2, new_df)
    
    return(new_df)
  }
  
  
  ### min-max based normalization
  normalization <- function(df) {
    
    maxV <- max(df)
    
    alpha <- 10 / maxV
    
    return(df * alpha)
    
  }
  
  
  ### function for circos
  makeCircos <- function(corMat, nodeType, fileName, outputPath) {
  
    ### set parameters for circos
    seg.num <- nrow(corMat)
    sample.num <- 10
    seg.name <- rownames(corMat)
    
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
    
    
    ### set db
    db <- segAnglePo(seg.f, seg=seg.name)
    
    ### set link.pg.v
    link.pg.v <- data.frame(matrix(0, length(which(corMat > 0)), 6))
    colnames(link.pg.v) <- c("seg1", "start1", "end1", "seg2", "start2", "end2")
    
    cnt <- 1
    for(i in 1:seg.num) {
      for(j in i:seg.num) {
        if(corMat[i,j] > 0) {
          link.pg.v[cnt,1] <- seg.name[i]
          link.pg.v[cnt,2] <- (sample.num/2)-(corMat[i,j]/2)
          link.pg.v[cnt,3] <- (sample.num/2)+(corMat[i,j]/2)
          link.pg.v[cnt,4] <- seg.name[j]
          link.pg.v[cnt,5] <- (sample.num/2)-(corMat[i,j]/2)
          link.pg.v[cnt,6] <- (sample.num/2)+(corMat[i,j]/2)
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
    circos(xc=450, yc=390, R=290, cir=db, W=40, mapping=link.pg.v, type="link.pg", lwd=2, col=sample(line_colors,nrow(link.pg.v), replace = TRUE))
    dev.off()
    
  }
  
  
  ### data1
  new_d1 <- makeSymmetricMat(d1)
  new_d1 <- normalization(new_d1)
  n <- c(rep("A", nrow(d1)), rep("B", ncol(d1)))
  names(n) <- rownames(new_d1)
  fileName <- "data1"
  makeCircos(new_d1, n, fileName, paste0(outputDir, fileName, ".pdf"))
  
  ### data2
  new_d2 <- makeSymmetricMat(d2)
  new_d2 <- normalization(new_d2)
  n <- c(rep("A", nrow(d2)), rep("B", ncol(d2)))
  names(n) <- rownames(new_d2)
  fileName <- "data2"
  makeCircos(new_d2, n, fileName, paste0(outputDir, fileName, ".pdf"))
  
  ### data3
  new_d3 <- makeSymmetricMat(d3)
  new_d3 <- normalization(new_d3)
  n <- c(rep("A", nrow(d3)), rep("B", ncol(d3)))
  names(n) <- rownames(new_d3)
  fileName <- "data3"
  makeCircos(new_d3, n, fileName, paste0(outputDir, fileName, ".pdf"))
  
}
