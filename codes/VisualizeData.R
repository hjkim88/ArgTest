###
#   File name : VisualizeData.R
#   Author    : Hyunjin Kim
#   Date      : Jan 2, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Performing TSNE to the datsets to see data status
#
#   Instruction
#               1. Source("VisualizeData.R")
#               2. Run the function "visualize()" - specify the input file (expression) directory and output path
#               3. t-SNE plots will be generated in the output path
#
#   Example
#               > source("The_directory_of_VisualizeData.R/VisualizeData.R")
#               > visualize(dataPath="./data/argtest.csv", outputPath="./results/Arg_tsne.pdf")
###


visualize <- function(dataPath="./data/argtest.csv", outputPath="./results/Arg_tsne.pdf") {
  
  ### load library
  if(!require(Rtsne)) {
    install.packages("Rtsne")
    library(Rtsne)
  }
  
  
  ### load dataset
  d <- read.csv(dataPath, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### keep type of the rows
  t <- d$Type
  names(t) <- rownames(d)
  d <- d[,-1]
  
  
  ### set sample and cluster information
  a <- cbind(as.character(colnames(d)), c("C1", "C1", "C2", "C2", "C2", "C2", "C3", "C3", "C3", "C3"))
  colnames(a) <- c("name", "cluster")
  a <- as.data.frame(a)
  
  
  ### run t-SNE
  set.seed(1234)
  tsne <- Rtsne(t(d), perplexity = 3)
  
  
  ### set color information
  colors = rainbow(length(unique(a$cluster)))
  names(colors) = unique(a$cluster)
  
  
  ### plot t-sne
  pdf(outputPath)
  plot(tsne$Y, col="white", xlab="tsne1", ylab="tsne2", main = "t-SNE plot", xlim = c(-150, 200))
  text(tsne$Y, labels=colnames(d), col=colors[a$cluster])
  legend("topright", legend = unique(a$cluster), col = colors[unique(a$cluster)], pch = 15)
  dev.off()
  
}
