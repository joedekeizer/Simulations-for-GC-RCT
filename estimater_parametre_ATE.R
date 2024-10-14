#################################################################
###### Get the empircal values
#################################################################
library(dplyr)
library(doParallel)

rm(list = ls())

gc()

path = "/home/joedekeizer/Simulations/"
Size.Effect = "log3"
sizeEffect = log(3)
N.effectif = 200

N.start = 1
N.stop = 1000000


setwd(path)

nb.cluster = parallel::detectCores() - 1 #- round(0.1 * parallel::detectCores())

if(.Platform[[1]] == "windows") {cl <- makeCluster(nb.cluster, type="PSOCK") } else{ cl <- makeCluster(nb.cluster, type="FORK") }
registerDoParallel(cl)
if(.Platform[[1]] == "windows") {clusterEvalQ(cl, {library(dplyr);library(MASS);library(splines);library(caret);library(SuperLearner);library(glmnet);library(cvAUC)})}



# library(doSNOW)
# cl <- makeSOCKcluster(90)
# registerDoSNOW(cl)




calcul.logOR = function(p0, p1)
{
  return(log((p1*(1-p0))/(p0*(1-p1))))
}
sim.base = function(N, beta1.t)
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
  
  
  return(data.obs)
}
res.sim = function(i, data.obs)
{
  p0 = sum(data.obs$outcome[data.obs$ttt == 0]/sum(1*I(data.obs$ttt == 0)))
  p1 = sum(data.obs$outcome[data.obs$ttt == 1]/sum(1*I(data.obs$ttt == 1)))
  
  # mod = summary(glm(outcome ~ ttt, family = binomial(link = "logit"), data = data.obs))
  # mod = prop.test(sum(data.obs$outcome[data.obs$ttt == 0]), sum(1*I(data.obs$ttt == 0)), correct = F)
  
  B = 1000
  p0.boot = p1.boot = logOR.boot = delta.boot = NA
  for(b in 1:B)
  {
    boot = data.obs[sample(1:N.effectif, size = N.effectif, replace = T),]
    
    p0.boot[b] = sum(boot$outcome[boot$ttt == 0]/sum(1*I(boot$ttt == 0)))
    p1.boot[b] = sum(boot$outcome[boot$ttt == 1]/sum(1*I(boot$ttt == 1)))
    logOR.boot[b] = calcul.logOR(p0 = p0.boot[b], p1 = p1.boot[b])
    delta.boot[b] = p1.boot[b] - p0.boot[b]
  }
  
  
  res = data.frame(
    iter = i,
    n0 = sum(1*I(data.obs$ttt == 0)),
    n1 = sum(1*I(data.obs$ttt == 1)),
    p0 = p0,
    p1 = p1,
    logOR = calcul.logOR(p0 = p0, p1 = p1),
    delta = p1 - p0,
    
    # sd.p0 = sqrt((p0*(1 - p0))/sum(1*I(data.obs$ttt == 0))),
    # sd.p1 = sqrt((p1*(1 - p1))/sum(1*I(data.obs$ttt == 1))),
    # sd.logOR = mod$coefficients[2,2],
    
    sd.p0 = sd(p0.boot),
    sd.p1 = sd(p1.boot),
    sd.logOR = sd(logOR.boot),
    sd.delta = sd(delta.boot)
  )
  
  return(res)
}


res <-  foreach(i = N.start:N.stop, .combine = rbind, .inorder = FALSE, .verbose = T) %dopar%
  {res.sim(i, data.obs = sim.base(N = N.effectif, beta1.t = sizeEffect))}


# registerDoSEQ()
stopCluster(cl)

save(res, file = paste0("res_n", N.effectif, "_ate", Size.Effect, ".Rdata"))

head(res)

print(mean(exp(res$logOR)))




