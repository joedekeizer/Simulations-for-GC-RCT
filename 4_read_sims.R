################################################################################
###### (4) Summarize the estimated results
################################################################################
### This code summarizes the results from the different methods on the 10000 simulated datasets into a table format (used to create Tables 1 and 2)
### It will save a .tex file (Latex) named Latex_complex_nXX_ateYY in the specified work directory where XX is the sample size (N.effectif) and YY the effect size (Size.Effect)
################################################################################
library(dplyr)
library(stringr)
library(purrr)
library(xtable)

### Initialization of the parameters
SizeEffect <- "log(3)" #The size effect (character)
N.effectif <- 200 #The sample size (numeric)
path <- "/home/Simulations/" #Work directory
path <- "C:\\Users\\Joe\\Documents\\GitHub\\Simulations-for-GC-RCT2\\" #Work directory
setwd(path)


################################################################################

### Functions used for calculating the different performance criteria (formulas detailed in the supplementary material)
calcul.esd <- function(mean.parameter, estimand)
{
  return(((mean.parameter - estimand)**2) %>% sum %>% `/`(length(mean.parameter)-1) %>% sqrt)
}
calcul.asd <- function(sd.parameter)
{
  return(mean(sd.parameter))
}
calcul.mean.bias <- function(parameter, estimand)
{
  return((parameter - estimand) %>% mean %>% `*`(100))
}
calcul.root.mean.square.error <- function(parameter, estimand)
{
  return(((parameter - estimand)**2) %>% mean %>% sqrt)
}
calcul.variance.estimation.bias <- function(mean.parameter, sd.parameter, estimand)
{
  return(((calcul.asd(sd.parameter) - calcul.esd(mean.parameter, estimand)) / calcul.esd(mean.parameter, estimand)) %>% `*`(100))
}
calcul.coverage <- function(borne.inf, borne.sup, estimand)
{
  return(mean(1*I(estimand > borne.inf & estimand < borne.sup)) %>% `*`(100))
}
calcul.puissance <- function(borne.inf, borne.sup)
{
  return(mean(1*I(borne.inf > 0) + 1*I(borne.sup < 0)) %>% `*`(100))
}
calcul.rss <- function(mean.parameter, sd.parameter, mean.parameter.unadjusted, sd.parameter.unadjusted)
{
  zref <- mean(mean.parameter.unadjusted)/mean(sd.parameter.unadjusted)
  zadj <- mean(mean.parameter)/mean(sd.parameter)
  return((1 - ((zref/zadj)**2)) %>% `*`(100))
}


### Load the theoretical estimates used to calculate bias for the different methods
Size.Effect <- gsub("[[:punct:]]", "", SizeEffect) #Removing any special character for folder name
load(paste0("res_n", N.effectif, "_ate", Size.Effect, ".Rdata"))
res_final <- res
res_final <- res_final[!(res_final$logOR=="-Inf" | res_final$logOR=="Inf"),] #In the very rare case where p0 = 0 or p1 = 1
rm(res)


### Results of the simulated data
pathresults = paste0(path, "/results_complex_n", N.effectif, "_ate", Size.Effect,"/")
liste = list.files(pathresults)


### Combining all the results together of each dataset
res = NULL
for(i in 1:length(liste))
{
  temp = read.delim(paste0(pathresults,liste[i]), header = T, sep = ";")
  res = rbind(res, temp)
}
rm(i, temp)


### Save the results as a table in a dataframe format
(dfl <- data.frame(
  Method = c("Unadjusted",
             "Elasticnet",
             "Lasso",
             "Neural network",
             "Support vector machine",
             "Super learner"
  ),
  
  ### Mean bias p0
  MBp0 = c(calcul.mean.bias(parameter = res$mean.p0.unadjusted.bcv, estimand = mean(res_final$p0)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p0.elasticnet.bcv, estimand = mean(res_final$p0)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p0.lasso.bcv, estimand = mean(res_final$p0)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p0.mlp.bcv, estimand = mean(res_final$p0)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p0.svmRadial.bcv, estimand = mean(res_final$p0)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p0.sl.bcv[which(abs(res$mean.p0.sl.bcv) < 100)], estimand = mean(res_final$p0)) %>% round(2)
  ),
  
  ### Mean bias p1
  MBp1 = c(calcul.mean.bias(parameter = res$mean.p1.unadjusted.bcv, estimand = mean(res_final$p1)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p1.elasticnet.bcv, estimand = mean(res_final$p1)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p1.lasso.bcv, estimand = mean(res_final$p1)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p1.mlp.bcv, estimand = mean(res_final$p1)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p1.svmRadial.bcv, estimand = mean(res_final$p1)) %>% round(2),
           calcul.mean.bias(parameter = res$mean.p1.sl.bcv[which(abs(res$mean.p1.sl.bcv) < 100)], estimand = mean(res_final$p1)) %>% round(2)
  ),
  
  ### Mean bias log mOR
  MBlogmOR = c((mean(res$mean.logOR.unadjusted.bcv) - mean(res_final$logOR)) %>% round(4),
             (mean(res$mean.logOR.elasticnet.bcv) - mean(res_final$logOR)) %>% round(4),
             (mean(res$mean.logOR.lasso.bcv) - mean(res_final$logOR)) %>% round(4),
             (mean(res$mean.logOR.mlp.bcv) - mean(res_final$logOR)) %>% round(4),
             (mean(res$mean.logOR.svmRadial.bcv) - mean(res_final$logOR)) %>% round(4),
             (mean(res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)]) - mean(res_final$logOR)) %>% round(4)
  ),
  
  ### Mean bias dekta
  MBdelta = c(calcul.mean.bias(parameter = res$mean.delta.unadjusted.bcv, estimand = mean(res_final$delta)) %>% round(2),
              calcul.mean.bias(parameter = res$mean.delta.elasticnet.bcv, estimand = mean(res_final$delta)) %>% round(2),
              calcul.mean.bias(parameter = res$mean.delta.lasso.bcv, estimand = mean(res_final$delta)) %>% round(2),
              calcul.mean.bias(parameter = res$mean.delta.mlp.bcv, estimand = mean(res_final$delta)) %>% round(2),
              calcul.mean.bias(parameter = res$mean.delta.svmRadial.bcv, estimand = mean(res_final$delta)) %>% round(2),
              calcul.mean.bias(parameter = res$mean.delta.sl.bcv[which(abs(res$mean.delta.sl.bcv) < 100)], estimand = mean(res_final$delta)) %>% round(2)
  ),
  
  ### Root mean square error (RMSE)
  RMSE = c(calcul.root.mean.square.error(parameter = res$mean.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
           calcul.root.mean.square.error(parameter = res$mean.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
           calcul.root.mean.square.error(parameter = res$mean.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
           calcul.root.mean.square.error(parameter = res$mean.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
           calcul.root.mean.square.error(parameter = res$mean.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
           calcul.root.mean.square.error(parameter = res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
  ),
  
  ### Variance estimation bias (VEB)
  VEB = c(calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.unadjusted.bcv, sd.parameter = res$sd.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
          calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.elasticnet.bcv, sd.parameter = res$sd.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
          calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.lasso.bcv, sd.parameter = res$sd.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
          calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.mlp.bcv, sd.parameter = res$sd.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
          calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.svmRadial.bcv, sd.parameter = res$sd.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
          calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], sd.parameter = res$sd.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
  ),
  
  ### Coverage
  Coverage = c(calcul.coverage(borne.inf = res$inf.logOR.unadjusted.bcv, borne.sup = res$sup.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.coverage(borne.inf = res$inf.logOR.elasticnet.bcv, borne.sup = res$sup.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.coverage(borne.inf = res$inf.logOR.lasso.bcv, borne.sup = res$sup.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.coverage(borne.inf = res$inf.logOR.mlp.bcv, borne.sup = res$sup.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.coverage(borne.inf = res$inf.logOR.svmRadial.bcv, borne.sup = res$sup.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.coverage(borne.inf = res$inf.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], borne.sup = res$sup.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
  ),
  
  ### Type I or II error according to the Size.Effect
  #If the size effect correspond to a mOR = 1 (under null hypothesis) then type I error is the calcul.puissance function
  #Else 100 - calcul.puissance function will return the type II error
  Error = c(ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.unadjusted.bcv, borne.sup = res$sup.logOR.unadjusted.bcv) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.unadjusted.bcv, borne.sup = res$sup.logOR.unadjusted.bcv) %>% round(2)),
            ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.elasticnet.bcv, borne.sup = res$sup.logOR.elasticnet.bcv) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.elasticnet.bcv, borne.sup = res$sup.logOR.elasticnet.bcv) %>% round(2)),
            ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.lasso.bcv, borne.sup = res$sup.logOR.lasso.bcv) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.lasso.bcv, borne.sup = res$sup.logOR.lasso.bcv) %>% round(2)),
            ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.mlp.bcv, borne.sup = res$sup.logOR.mlp.bcv) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.mlp.bcv, borne.sup = res$sup.logOR.mlp.bcv) %>% round(2)),
            ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.svmRadial.bcv, borne.sup = res$sup.logOR.svmRadial.bcv) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.svmRadial.bcv, borne.sup = res$sup.logOR.svmRadial.bcv) %>% round(2)),
            ifelse(Size.Effect == "log09729" | Size.Effect == "log1", calcul.puissance(borne.inf = res$inf.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], borne.sup = res$sup.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)]) %>% round(2),
                   100 - calcul.puissance(borne.inf = res$inf.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], borne.sup = res$sup.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)]) %>% round(2))
  ),
  
  ### Reduction of sample size (RSS)
  RSS = c(0,
          calcul.rss(mean.parameter = res$mean.logOR.elasticnet.bcv, sd.parameter = res$sd.logOR.elasticnet.bcv,
                     mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2),
          calcul.rss(mean.parameter = res$mean.logOR.lasso.bcv, sd.parameter = res$sd.logOR.lasso.bcv,
                     mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2),
          calcul.rss(mean.parameter = res$mean.logOR.mlp.bcv, sd.parameter = res$sd.logOR.mlp.bcv,
                     mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2),
          calcul.rss(mean.parameter = res$mean.logOR.svmRadial.bcv, sd.parameter = res$sd.logOR.svmRadial.bcv,
                     mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2),
          calcul.rss(mean.parameter = res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], sd.parameter = res$sd.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)],
                     mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2)
  )
))


### Formatting the number of decimals (keep leading zeros)
dfl[,2] <- paste0(sprintf("%.2f", dfl[,2]), "%")
dfl[,3] <- paste0(sprintf("%.2f", dfl[,3]), "%")
dfl[,4] <- sprintf("%.4f", dfl[,4])
dfl[,5] <- paste0(sprintf("%.2f", dfl[,5]), "%")
dfl[,6] <- formatC(dfl[,6], format = "f", digits = 2)
dfl[,7] <- paste0(sprintf("%.2f", dfl[,7]), "%")
dfl[,8] <- paste0(sprintf("%.2f", dfl[,8]), "%")
dfl[,9] <- paste0(sprintf("%.2f", dfl[,9]), "%")
dfl[,10] <- paste0(sprintf("%.2f", dfl[,10]), "%")


################################################################################

### Table with all the results for the specific scenario
print(dfl)

### Save the result table in a .tex format
print(xtable(dfl, type = "latex"), file = paste0(paste0("Latex_complex_n", N.effectif, "_ate", Size.Effect,".tex")),include.rownames=FALSE)
