###
#   File name : MakeCircularPlot.R
#   Author    : Hyunjin Kim
#   Date      : Jan 3, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate Spearman corrleation and draw a circular plot based on the correlation
#
#   Instruction
#               1. Source("MakeCircularPlot.R")
#               2. Run the function "makeCirPlot" - specify the input files (microbe & antibiotic) and output directory
#               3. The circular plot will be generated in the output directory
#
#   Example
#               > source("The_directory_of_MakeCircularPlot.R/MakeCircularPlot.R")
#               > makeCirPlot(dataPath="./data/argtestre_re.csv", outputDir="./results/")
###


makeCirPlot <- function(dataPath="./data/argtestre_re.csv", outputDir="./results/") {
  
  ### load CorNetCirPlot.R
  source("./codes/CorNetCirPlot.R")
  
  
  ### load datasets
  d <- read.csv(dataPath, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  
  ### keep type of the rows
  t <- d$Type
  names(t) <- rownames(d)
  d <- d[,-1]
  
  
  ### HRSD, CA, DEMON & EBBR
  d_hrsd <- d[,1:4]
  d_ca <- d[,9:10]
  d_de <- d[,5:8]
  
  
  ### cirplot
  cirPlot(d, t, outputPath = paste0(outputDir, "CorNet_cirplot.pdf"))
  cirPlot(d_hrsd, t, outputPath = paste0(outputDir, "CorNet_cirplot_HRSD.pdf"))
  cirPlot(d_ca, t, outputPath = paste0(outputDir, "CorNet_cirplot_CA.pdf"))
  cirPlot(d_de, t, outputPath = paste0(outputDir, "CorNet_cirplot_DE.pdf"))
  
}

