# tcontam and tcontamL type 1 error analysis

library(decontam)
library(tidyverse)
library(dirmult)
library(foreach)
library(doParallel)

load("oral_data_new.Rdata")
out <- as.matrix(otu)
otu <- out[, order(colSums(otu), decreasing = T)]

source("decontam_tcontamL.R")

# set the biomass
total1 <- runif(10, min = 95, max = 120)

# estimate the DM parameters
set.seed(123)
fit <- dirmult::dirmult(otu)

pi     <- fit$p
theta  <- fit$theta
theta1 <- theta 


otuname <- names(pi)
topotuname1 <- names(sort(pi, decreasing = T)[1:200])

toppi1 <- pi[topotuname1] /  sum(pi[topotuname1])
shape1 <- toppi1 * ((1 - theta1) /  theta1)

myiter <- 2000


registerDoParallel(detectCores()-1)


alltype1 <- foreach(i = 1:myiter, .inorder = FALSE)%dopar%{
  
  print(i)
  set.seed(i*100)
  prob1 <- rdirichlet(10, shape1)

  
  biomass1 <- sweep(prob1, 1, total1, FUN = "*")

  mytable <- cbind(biomass1)
  mytable <- mytable/rowSums(mytable)
  colnames(mytable) <- paste("Seq", 1:ncol(mytable))
  
  conc <- total1 
  
  
  ocf_original <- isContaminant(mytable, conc=conc, method = "frequency")
  ocf_test <- isContaminant(mytable, conc = conc, method = "frequency", frequency.method = "test")
  ocf_test_low_0.05 <- isContaminant(mytable, conc = conc, method = "frequency", frequency.method = "test", lowbiomass = T, lowthrehold = 0.05)
  
  alliter <- data.frame(seqID = colnames(mytable),
                        prev = ocf_original$prev,
                        decontam = ocf_original$p,
                        tcontam = ocf_test$p,
                        tcontamL_0.05 = ocf_test_low_0.05$p)
  
  return(alliter)
  
}



save(alltype1, file = "simu_real_type1.Rdata")

dev.off()
