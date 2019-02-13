#calculate values of performance metrics 


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

setwd("C:/Users/anush/OneDrive/Research Lab/")
# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)  # helper function from STIR


read.data <- function(sim.type, n.samples, n.attributes, pct.signals, label, bias, training, holdout, validation, 
                      verbose){
  sim.type <- sim.type
  
  #read analysis file and get confusion matrix data 
  pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                       "pct.signals", pct.signals,
                       "num.attr", n.attributes, "num.samp", n.samples, sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="")
  dat <- read.csv(pec_simFile)

}

data <-read.data(sim.type = "mainEffect", n.samples = 100,
              n.attributes = 1000,
              pct.signals = 0.1, label ="class", bias = 0.8,
              training = 1/3,
              holdout = 1/3,
              validation = 1/3,
              verbose = FALSE)

