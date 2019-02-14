# simulate and write data 
# oder of execution 1

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


# INPUT VALUES 
letsSimulate <- T   # F to use previously simulated data
writeData <- T  # usually the same as letsSimulate
writeResults <- F

simulate.data <- function(sim.type, n.samples, n.attributes, pct.signals, label, bias, training, holdout, validation, 
                          verbose, replication){
  sim.type <- sim.type
  metric<- "manhattan"
  for (i in 1:replication){
    sim.data<- createSimulation(num.samples = n.samples,
                                num.variables = n.attributes,# attributes
                                pct.signals = pct.signals,# percentage of simulated values
                                label = label,
                                bias = bias,# strength of main and interaction effect in simulated values
                                pct.train = training,
                                pct.holdout = holdout,
                                pct.validation = validation,
                                sim.type = sim.type,
                                verbose = verbose)
    dat <- rbind(sim.data$train, sim.data$holdout)
    # if (i == 1 ){
    #   result.data <- data 
    # }else {
    #   rbind(result.data, dat)
    # }
    
    #simulate and save the data 
    pec_simFile <- paste(i,"pec_simulated", sim.type, "bias", bias, 
                         "pct.signals", pct.signals,
                         "num.attr", n.attributes, "num.samp", n.samples, sep = "_")
    pec_simFile <- paste(pec_simFile,".csv",sep="")
    #dat <- stack(as.data.frame(result.data))
    setwd("C:/Users/anush/OneDrive/Research Lab/Project/In Silico Feature Selection/data")
    write.csv(dat, file=pec_simFile, row.names = FALSE)
  }
}

simulate.data(sim.type = "mainEffect", n.samples = 100,
                               n.attributes = 1000,
                               pct.signals = 0.1, label ="class", bias = 0.8,
                               training = 1/3,
                               holdout = 1/3,
                               validation = 1/3,
                               verbose = FALSE, replication = 2)

# simulate.data(sim.type = "interactionErdos", n.samples = 100,
#               n.attributes = 1000,
#               pct.signals = 0.1, label ="class", bias = 0.4,
#               training = 1/2,
#               holdout = 1/2,
#               validation = 0,
#               verbose = FALSE, replication = 1)


  