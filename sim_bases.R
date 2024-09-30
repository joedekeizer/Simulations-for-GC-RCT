###### Simulation bases 


library(dplyr)
library(doSNOW)

rm(list = ls())

gc()

path = "/home/joedekeizer/Simulations/"
path = "C:/Users/Joe/Documents/Thesis/PRGM/Simulations/"
Size.Effect = "log3"
sizeEffect = log(3)
N.effectif = 350

parallel::detectCores()
cl <- makeSOCKcluster(15) 
registerDoSNOW(cl)

# 10k bases to simulate and save as Rdata (data.obs)
N.start = 1
N.stop = 10000

pathbases = paste0(path, "/complex_n", N.effectif, "_ate", Size.Effect, "/")

#Create the folder if it does not already exist
ifelse(!dir.exists(file.path(pathbases)), dir.create(file.path(pathbases)), FALSE)

setwd(pathbases)


sim.base = function(N, beta1.t, iter)
{
  beta0.x = -0.4
  beta1.x = log(2)
  
  beta0.o = -2
  beta1.o = log(2)
  
  
  .x1 <- rnorm(N, 0, 1) #1
  .x2 <- rnorm(N, beta0.x + beta1.x * .x1, 1) #2
  .x3 <- rnorm(N, beta0.x - beta1.x * .x1 - beta1.x * .x2, 1) #3
  .x5 <- rnorm(N, 0, 1) #4
  .x6 <- 1 * (rnorm(N, 0, 1) > 0.66) # prevalence ~25% #5
  .x7 <- 1 * (rnorm(N, beta0.x - beta1.x * .x5, 1) > (-0.40)) # prevalence ~50% #6
  .x8 <- rnorm(N, beta0.x - beta1.x * .x6, 1) #7
  .x9 <- 1 * (rnorm(N, beta0.x + beta1.x * .x7, 1) > (-0.80)) # prevalence ~75% (77%) #8
  .x10 <- rnorm(N, beta0.x + beta1.x * .x8, 1) #9
  .x11 <- rnorm(N, 0, 1) #10
  .x12 <- 1 * (rnorm(N, beta0.x + beta1.x * .x9, 1) > (0.84)) # prevalence ~25% #11
  .x14 <- rnorm(N, beta0.x - beta1.x * .x12 - beta1.x * .x11, 1) #12
  .x15 <- rnorm(N, beta0.x - beta1.x * .x12, 1) #13
  .x18 <- rnorm(N, 0, 1) #14
  .x19 <- 1 * (rnorm(N, 0, 1) > qnorm(0.75))  # prevalence ~25% #15
  .x20 <- 1 * (rnorm(N, 0, 1) > (0.66))  # prevalence ~25% #16
  .x21 <- rnorm(N, 0, 1) #17
  
  data.obs <- data.frame(x1 = .x1,
                         x2 = .x2,
                         x3 = .x3,
                         x5 = .x5,
                         x6 = .x6,
                         x7 = .x7,
                         x8 = .x8,
                         x9 = .x9,
                         x10 = .x10,
                         x11 = .x11,
                         x12 = .x12,
                         x14 = .x14,
                         x15 = .x15,
                         x18 = .x18,
                         x19 = .x19,
                         x20 = .x20,
                         x21 = .x21)
  data.obs$ttt = rbinom(N, size = 1, prob = 0.5)
  
  # mean(data.obs$ttt)
  
  bx <- beta0.o +
    beta1.o * (data.obs$x2 > -0.44) -
    beta1.o * data.obs$x3 + (beta1.o / 2) * (data.obs$x3^2) +
    beta1.o * data.obs$x6 +
    beta1.o * data.obs$x7 +
    beta1.o * data.obs$x10 +
    beta1.o * 0.5 * (data.obs$x11^2) -
    beta1.o * data.obs$x14 -
    beta1.o * (data.obs$x15 > -0.55) +
    beta1.o * data.obs$x18 +
    beta1.o * data.obs$x19 +
    beta1.o * 0.5 * data.obs$ttt*data.obs$x18 +
    beta1.t * data.obs$ttt
  
  
  
  pr.o.obs = plogis(bx)
  mean(pr.o.obs)
  
  data.obs$outcome <- rbinom(N, 1, prob = pr.o.obs)
  
  save(data.obs, file = paste0(pathbases,"sim_n", N.effectif, "_ate", Size.Effect,"_",formatC(iter, width = 5, format = "d", flag = "0"), ".Rdata"))
  
  return(data.obs)
}



foreach(i = N.start:N.stop,  .inorder = FALSE, .verbose = T) %dopar%
  {sim.base(N = N.effectif, beta1.t = sizeEffect, iter = i)}



# registerDoSEQ()
stopCluster(cl)
