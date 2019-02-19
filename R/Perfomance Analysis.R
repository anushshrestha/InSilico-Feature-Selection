# filter all result from analysis with top 10 < threshold value 
# calculate values of performance metrics 
# order of execution 3

knitr::opts_chunk$set(echo = TRUE, warning=FALSE)
rm(list = ls())

if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}

library(devtools)

if (!("privateEC" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/privateEC") # build_vignettes = TRUE)
}
if (!("stir" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/stir", build_vignettes = TRUE)
}
library(privateEC)  # used to simulate data
library(stir)

# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)  # helper function from STIR

build.read.filename <- function(i,sim.type, bias, pct.signals, n.attributes, n.samples){
  
  pec_simFile <- paste(i, "pec_simulated", sim.type, "bias", bias, 
                       "pct.signals", pct.signals,
                       "num.attr", n.attributes, "num.samp", n.samples, sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="") 
}

build.write.filename <- function(info.string){
  pec_simFile <- paste(info.string, sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="")
}

calculate.perform.metrices <- function(n.samples, n.attributes, threshold){
  currentPath <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(currentPath)
  setwd("../output")
  file.list <- list.files(getwd(), pattern="^[0-9]_[0-9a-zA-Z_. ]+.csv$")
  # t_sorted_relieff <- read.csv(file.list[2])
  # positive <- sum(t_sorted_relieff[,3]< 0.05)
  # negative <- 1000 - positive
  
  result.matrix <- matrix(0,nrow=length(file.list), ncol=5, byrow = TRUE)
  
  for(i in 1:length(file.list)){
    # while reading file row names can be taken from some fixed column 
    t_sorted_relieff <- read.csv(file.list[i], row.names = 1)
    # build confusion matrix
    positive <- sum(t_sorted_relieff[,3]< threshold)
    negative <- n.attributes - positive

    true.positive <- grep("^simvar[0-9]",rownames(t_sorted_relieff[1:positive,]), fixed = FALSE, value = TRUE)

    n.true.positive <- length(true.positive)
    false.positive <- grep("^var[0-9]", rownames(t_sorted_relieff[1:positive,]), fixed = FALSE, value = TRUE)
    n.false.positive <- length(false.positive)
    
    true.negative <-  grep("^var[0-9]",rownames(t_sorted_relieff[negative:n.samples,]), fixed = FALSE, value = TRUE)
    n.true.negative <- length(true.negative)
    false.negative <- grep("^simvar[0-9]",rownames(t_sorted_relieff[negative:n.samples,]), fixed = FALSE, value = TRUE)
    n.false.negative <- length(false.negative)
      
    # performance comparison
    tnr <- n.true.negative/(n.true.negative + n.false.positive)
    precision <- n.true.positive /(n.true.positive + n.false.positive)
    recall <- n.true.positive/(n.true.positive+ n.false.negative)
    
    # grab sim.type and analysis.type from file name 
    split.filename <-strsplit(file.list[i],split='_', fixed=TRUE) 
    analysis.type = split.filename[[1]][3]
    simulation.type = split.filename[[1]][6]

    result.matrix[i,] <- c(simulation.type, analysis.type, tnr, precision, recall)
    print(result.matrix[i,])
  }
  
  colnames(result.matrix) <- c("sim.type","analy.type","tnr","precision","recall")
  (result.matrix)
  dat <- as.data.frame(result.matrix)
  
  write.filename <- build.write.filename("performance metrics")
  
  currentPath <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(currentPath)
  setwd("../output")
  write.csv(dat, file = write.filename, row.names = FALSE) 
}

calculate.perform.metrices(n.samples = 100, n.attributes = 1000, threshold = 0.05)


