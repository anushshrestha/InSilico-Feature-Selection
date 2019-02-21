# analysis of data, calculate all the output of the analysis
# order of execution 2

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
install_github("insilico/npdr")
library(npdr)
library(glmnet)

# load other helper packages
packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
check.packages(packages)  # helper function from STIR

build.read.filename <- function(i,sim.type, bias, pct.signals, n.attributes, n.samples){
  
  pec_simFile <- paste(i, "pec_simulated", sim.type, "bias", bias, 
                       "pct.signals", pct.signals,
                       "num.attr", n.attributes, "num.samp", n.samples, sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="") 
}

build.write.filename <- function(i,analysis.type, sim.type, bias, pct.signals, n.attributes, n.samples){
  pec_simFile <- paste(i,"analyzed.with", analysis.type, "pec_simulated", sim.type, "bias", bias, 
                       "pct.signals", pct.signals,
                       "num.attr", n.attributes, "num.samp", n.samples,"analysis", sep = "_")
  pec_simFile <- paste(pec_simFile,".csv",sep="")
}

analyze.run <- function(data, analysis.type){
  if(analysis.type == "STIR SURF") {
    result <- STIR.SURF.analysis(data)
    }
  else if(analysis.type == "STIR MULTISURF") {
    result <- STIR.MULTISURF.analysis(data)
  }
  return (result)
}

# replication is equal to number of replication from simulate data 
analyze.data <- function(sim.type, n.samples, n.attributes, pct.signals, label, bias, training, holdout, validation, 
                          verbose, replication, analyze.todo){
  for (i in 1:replication){
    sim.type <- sim.type
    # path update
    read.filename <-build.read.filename(i,sim.type, bias, pct.signals, n.attributes, n.samples)
    currentPath <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(currentPath)
    setwd("../data")
    read.filename <-build.read.filename(i,sim.type, bias, pct.signals, n.attributes, n.samples)
    dat <- read.csv(read.filename)
    # for a given data, run multiple analysis as given 
    # after reading, prepare env for writing 
    currentPath <- dirname(rstudioapi::getSourceEditorContext()$path)
    #path update
    setwd(currentPath)
    setwd("../output")
    for (analyze.type in analyze.todo){
      t_sorted_surf <- analyze.run(data=dat, analysis.type= analyze.type)
      analyze.run(data = dat, analysis.type = analyze.type)
      write.filename <- build.write.filename(i,analyze.type,sim.type, bias, pct.signals, n.attributes, n.samples)
      write.csv(t_sorted_surf, file = write.filename, row.names = TRUE) 
    }
  }
}

# pass dat (single file) 
# main effect has i=1 and interaction has i=i+1
KNN.analysis <- function(dat){
  # MAIN EFFECT KNN 
  t_sorted_relieff <- list()
  RF.method = "relieff"
  metric <- "manhattan"
  k <- 16  # k=m/6 should be similar to MultiSURF
  i <- 1  # if you want to use k for loop

  class.label <- "class" 
  predictors.mat <- dat[, - which(colnames(dat) == class.label)]
  dat[, class.label] <- as.factor(dat[, class.label]) 
  pheno.class <- dat[, class.label]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)

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

STIR.SURF.analysis <- function(dat){
  # MAIN EFFECT SURF K = 0
  
  class.label <- "class" 
  predictors.mat <- dat[, - which(colnames(dat) == class.label)]
  dat[, class.label] <- as.factor(dat[, class.label]) 
  pheno.class <- dat[, class.label]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)
  
  t_sorted_relieff <- list()
  metric <- "manhattan"
  RF.method = "surf"
  k <- 0
  # k=m/6 should be similar to MultiSURF
  #i <- 1  # if you want to use k for loop
  neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.class, k = k, method = RF.method)
  results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
  t_sorted_surf <- results.list$STIR_T[, -3]
  colnames(t_sorted_surf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), sep=".")
  return(t_sorted_surf)
}

STIR.MULTISURF.analysis <- function(dat){
  
  class.label <- "class" 
  predictors.mat <- dat[, - which(colnames(dat) == class.label)]
  dat[, class.label] <- as.factor(dat[, class.label]) 
  pheno.class <- dat[, class.label]
  attr.names <- colnames(predictors.mat)
  num.samp <- nrow(dat)
  
  t_sorted_relieff <- list()
  RF.method = "multisurf"
  metric <- "manhattan"
  # let k=0 because multisurf does not use k
  neighbor.idx.observed<- find.neighbors(predictors.mat, pheno.class, k = 0, method = RF.method)
  results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = metric, method = RF.method)
  t_sorted_multisurf <- results.list$STIR_T[, -3]  # remove cohen-d
  colnames(t_sorted_multisurf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), sep=".")
  return(t_sorted_multisurf)
}

GLMSTIR.SURF.analysis <- function(dat){
  case.control.data <- dat
  n.samples.case.control <- dim(case.control.data)[1]
  pheno.case.control <- as.factor(case.control.data[,"class"])
  
  npdr.cc.results <- npdr("class", case.control.data, regression.type = "binomial", 
                          attr.diff.type = "numeric-abs", nbd.methods = "surf", 
                          nbd.metric = "manhattan", msurf.sd.frac = .5, 
                          padj.method = "bonferroni", verbose = T)
  colnames(npdr.cc.results) <- paste(c("pval.adj", "pval.att", "beta.raw.att",
                                       "beta.Z.att", "beta.0", "pval.0"), "glmstir", 
                                     sep = ".")
  return (npdr.cc.results)
}

GLMSTIR.MULTISURF.analysis <- function(dat){
  case.control.data <- dat
  n.samples.case.control <- dim(case.control.data)[1]
  pheno.case.control <- as.factor(case.control.data[,"class"])
  
  npdr.cc.results <- npdr("class", case.control.data, regression.type = "binomial", 
                          attr.diff.type = "numeric-abs", nbd.methods = "multisurf", 
                          nbd.metric = "manhattan", msurf.sd.frac = .5, 
                          padj.method = "bonferroni", verbose = T)
  colnames(npdr.cc.results) <- paste(c("pval.adj", "pval.att", "beta.raw.att",
                                       "beta.Z.att", "beta.0", "pval.0"), "glmstir", 
                                     sep = ".")
  return (npdr.cc.results)
}

univariate.analysis <- function(dat){
  case.control.data <- dat
  n.samples.case.control <- dim(case.control.data)[1]
  pheno.case.control <- as.factor(case.control.data[,"class"])
  univariate.cc.results <- uniReg(outcome="class", dataset=case.control.data, regression.type="binomial")
  return(univariate.cc.results)
}

glmnet.analysis <- function(dat){
  case.control.data <- dat
  n.samples.case.control <- dim(case.control.data)[1]
  pheno.case.control <- as.factor(case.control.data[,"class"])
  
  
}


# reading data from file 
analyzed.data <- analyze.data(sim.type = "mainEffect", n.samples = 100, 
                              n.attributes = 1000, 
                              pct.signals = 0.1, label ="class", bias = 0.8, 
                              training = 1/3,
                              holdout = 1/3,
                              validation = 1/3, 
                              verbose = FALSE,
                              replication = 2,
                              analyze.todo=c("STIR SURF",
                                             "STIR MULTISURF",
                                             "GLMSTIR SURF",
                                             "GLMSTIR MULTISURF",
                                             "univariate",
                                             "glmnet"
                                             )
                              )

# analyzed.data <- analyze.data(sim.type = "interactionErdos", n.samples = 100,
#               n.attributes = 1000,
#               pct.signals = 0.1, label ="class", bias = 0.4,
#               training = 1/2,
#               holdout = 1/2,
#               validation = 0,
#               verbose = FALSE, replication = 1)
