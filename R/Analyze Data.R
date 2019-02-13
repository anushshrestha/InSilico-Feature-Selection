#analysis of data, calculate values of confusion matrix

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
  
  pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                       "pct.signals", pct.signals,
                       "num.attr", n.attributes, "num.samp", n.samples, sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="")
  dat <- read.csv(pec_simFile)
}

# reading data from file 
all.dat <- read.data(sim.type = "mainEffect", n.samples = 100, 
              n.attributes = 1000, 
              pct.signals = 0.1, label ="class", bias = 0.8, 
              training = 1/3,
              holdout = 1/3,
              validation = 1/3, 
              verbose = FALSE)




# run loop and get single data 
for (j in 1:nrow(all.dat)){
  # MAIN EFFECT KNN 
  t_sorted_relieff <- list()
  RF.method = "relieff"
  metric <- "manhattan"
  k <- 16  # k=m/6 should be similar to MultiSURF
  i <- 1  # if you want to use k for loop
  
  dat <- all.dat[1,]
  class.label <- "class" 
  predictors.mat <- dat[, - which(colnames(dat) == class.label)]
  dat[, class.label] <- as.factor(dat[, class.label]) 
   <- dat[, class.label]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)
  
  browser()
  neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = k, method = RF.method)
  results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
  t_sorted_relieff[[1]] <- results.list$STIR_T[, -3]
  colnames(t_sorted_relieff[[i]]) <- paste(c("t.stat", "t.pval", "t.pval.adj"), k, sep=".")
  (t_sorted_relieff[[1]][1:10,])
  #build confusion matrix
  positive <- sum(t_sorted_relieff[[1]][,3]< 0.05)
  negative <- n.attributes - positive
  
  true.positive <- grep("^simvar[0-9]",rownames(t_sorted_relieff[[i]][1:positive,]), fixed = FALSE, value = TRUE)
  n.true.positive <- length(true.positive)
  false.positive <- grep("^var[0-9]", rownames(t_sorted_relieff[[i]][1:positive,]), fixed = FALSE, value = TRUE)
  n.false.positive <- length(false.positive)
  
  true.negative <-  grep("^var[0-9]",rownames(t_sorted_relieff[[i]][negative:n.samples,]), fixed = FALSE, value = TRUE)
  n.true.negative <- length(true.negative)
  false.negative <- grep("^simvar[0-9]",rownames(t_sorted_relieff[[i]][negative:n.samples,]), fixed = FALSE, value = TRUE)
  n.false.negative <- length(false.negative)
  
  result.data <- rbind(result.data, c(n.true.positive,n.false.positive,n.true.negative,n.false.negative))
}

pec_simFile <- paste("pec_simulated", sim.type, "bias", bias, 
                     "pct.signals", pct.signals,
                     "num.attr", n.attributes, "num.samp", n.samples,"analysis", sep = "_")
pec_simFile <- paste(pec_simFile,".csv",sep="")
dat <- stack(as.data.frame(result.data))
write.csv(dat, file=pec_simFile)