library(decontam)
library(tidyverse)
library(dirmult)
library(foreach)
library(doParallel)

load("oral_data_new.Rdata")
out <- as.matrix(otu)
otu <- out[, order(colSums(otu), decreasing = T)]


source("decontam_tcontamL.R")

set.seed(123)
total1 <- runif(10, min = 5, max = 40)
total2 <- rep(5, 10)

set.seed(123)
fit <- dirmult::dirmult(otu)

pi <- fit$p
theta <- fit$theta



otuname <- names(pi)
topotuname1 <- names(sort(pi, decreasing = T)[1:200])
topotuname2 <- names(sort(pi, decreasing = T)[201:250])

toppi1 <- pi[topotuname1] /  sum(pi[topotuname1])
toppi2 <- pi[topotuname2] /  sum(pi[topotuname2])




theta1 <- theta
shape1 <- toppi1 * ((1 - theta1) /  theta1)

theta2 <- seq(0.001, 0.004, 0.0005)


shape2 <- lapply(theta2, function(x){toppi2 * ((1 - x) /  x)})



myiter <- 500


registerDoParallel(detectCores()-1)


allpower_verylow_all <- list()
for(jj in 1:length(shape2)){
  
  allpower_low <- foreach(i = 1:myiter, .inorder = FALSE)%dopar%{
    
    
    print(i)
    set.seed(i*100)
    prob1 <- rdirichlet(10, shape1)
    prob2 <- rdirichlet(10, shape2[[jj]])
    
    
    biomass1 <- sweep(prob1, 1, total1, FUN = "*")
    biomass2 <- sweep(prob2, 1, total2, FUN = "*")
    
    mytable <- cbind(biomass1, biomass2)
    mytable <- mytable/rowSums(mytable)
    colnames(mytable) <- paste("Seq", 1:ncol(mytable))
    
    conc <- total1 + total2
    
    
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
  
  allpower_verylow_all[[jj]] <- allpower_low
  
}



save(allpower_verylow_all, file = "simu_real_power_verylow_all_theta.Rdata")


