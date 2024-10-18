#################################################################
###### (1) Simulations to create the data sets
#################################################################
### This code will create a folder in the specified work directory named /complex_nXX_ateYY where XX is the sample size (N.effectif) and YY the effect size (Size.Effect)
### In the folder 10000 simulated data sets will be created and saved as .Rdata, according to the relations described in the sim.base function
#################################################################
library(dplyr)
library(doParallel)

### Initialization
path <- "/home/Simulations/" #Work directory
Size.Effect <- "log3" #The beta4 coefficient from Supplementary Table A1 for folder names (character format)
sizeEffect <- log(3) #Same beta4 coefficient but numeric for the simulations (numeric)
N.effectif <- 200 #Sample size (numeric)
N.stop <- 10000 #10k data sets to simulate and save as Rdata (data.obs)
setwd(path)

### Parallelisation of the code (PSOCK if Windows and FORK if Linux environment)
nb.cluster <- parallel::detectCores() - 1 #All cores minus one
if(.Platform[[1]] == "windows") {cl <- makeCluster(nb.cluster, type="PSOCK") } else{ cl <- makeCluster(nb.cluster, type="FORK") }
registerDoParallel(cl)
if(.Platform[[1]] == "windows") {clusterEvalQ(cl, {library(dplyr);library(MASS);library(splines);library(caret);library(SuperLearner);library(glmnet);library(cvAUC)})}


### Create the folder for the data if it does not already exist
pathbases <- paste0(path, "/complex_n", N.effectif, "_ate", Size.Effect, "/")
ifelse(!dir.exists(file.path(pathbases)), dir.create(file.path(pathbases)), FALSE)


#################################################################

### Function to simulate the data
sim.base <- function(N, beta1.t, iter)
{
  ### Coefficient from Supplementary Table A1
  beta0.x <- -0.4
  beta1.x <- log(2)
  
  beta0.o <- -2
  beta1.o <- log(2)
  
  ### Simulated variables from the complex scenario Figure 1
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
  
  ### Causal relations between the variables and the outcome from the complex scenario Figure 1
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
  
  pr.o.obs <- plogis(bx)
  
  data.obs$outcome <- rbinom(N, 1, prob = pr.o.obs)
  
  ### Save each data set as data.obs in a .Rdata file
  save(data.obs, file = paste0(pathbases,"sim_n", N.effectif, "_ate", Size.Effect,"_",formatC(iter, width = 5, format = "d", flag = "0"), ".Rdata"))
  
  return(data.obs)
}


#################################################################

### Running the function with parallel processing (not necessary but recommended)
if(.Platform[[1]] == "windows") {
  ### Windows
  foreach(i = 1:N.stop,  .inorder = FALSE, .verbose = T,
          .export=ls(envir=globalenv()),
          .packages = c("dplyr")
  ) %dopar% {sim.base(N = N.effectif, beta1.t = sizeEffect, iter = i)}
} else {
  ### Linux
  foreach(i = 1:N.stop,  .inorder = FALSE, .verbose = T
  ) %dopar% {sim.base(N = N.effectif, beta1.t = sizeEffect, iter = i)}
}

stopCluster(cl)