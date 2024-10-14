library(dplyr)
library(stringr)
library(purrr)
library(xtable)


rm(list = ls())

calcul.esd = function(mean.parameter, estimand)
{
  return(((mean.parameter - estimand)**2) %>% sum %>% `/`(length(mean.parameter)-1) %>% sqrt)
}
calcul.asd = function(sd.parameter)
{
  return(mean(sd.parameter))
}
calcul.mean.bias = function(parameter, estimand)
{
  return((parameter - estimand) %>% mean %>% `*`(100))
}
calcul.root.mean.square.error = function(parameter, estimand)
{
  return(((parameter - estimand)**2) %>% mean %>% sqrt)
}
calcul.variance.estimation.bias = function(mean.parameter, sd.parameter, estimand)
{
  return(((calcul.asd(sd.parameter) - calcul.esd(mean.parameter, estimand)) / calcul.esd(mean.parameter, estimand)) %>% `*`(100))
}
calcul.coverage = function(borne.inf, borne.sup, estimand)
{
  return(mean(1*I(estimand > borne.inf & estimand < borne.sup)) %>% `*`(100))
}
calcul.puissance = function(borne.inf, borne.sup)
{
  return(mean(1*I(borne.inf > 0) + 1*I(borne.sup < 0)) %>% `*`(100))
}
calcul.rss = function(mean.parameter, sd.parameter, mean.parameter.unadjusted, sd.parameter.unadjusted)
{
  zref = mean(mean.parameter.unadjusted)/mean(sd.parameter.unadjusted)
  zadj = mean(mean.parameter)/mean(sd.parameter)
  return((1 - ((zref/zadj)**2)) %>% `*`(100))
}




path = "C:/Users/Joe/Documents/Thesis/PRGM/Simulations/"
setwd(path)



# Size.Effect = "log09729"
# Size.Effect = "log1.5"
# Size.Effect = "log3"
# N.effectif = 60
# N.effectif = 100
# N.effectif = 200




SE = c("log09729","log1.5","log3")
NE = c(60,100,200)

SE = c("log3")
NE = c(350)
for (s in SE) {
  for (j in NE) {
    
    Size.Effect = s
    N.effectif = j
    
    
    load(paste0("res_n", N.effectif, "_ate", Size.Effect, ".Rdata"))
    res_final = res
    res_final = res_final[!(res_final$logOR=="-Inf" | res_final$logOR=="Inf"),] #CFHERE need to recheck why many NAN for sd.logOR
    rm(res)
    
    pathsims=paste0(path, "/resultats_complex_n", N.effectif, "_ate", Size.Effect,"/")
    liste = list.files(pathsims)
    
    
    #Only using the first 10000 simulation results (11000 simulated results due to low n having no events in some bootstraps)
    res = NULL
    
    for(i in 1:length(liste))
    {
      temp = read.delim(paste0(pathsims,liste[i]), header = T, sep = ";")
      
      #Removing variables from results that are not used in the tables (new simple code does not contain these variables):
      temp=temp[,!grepl("app",names(temp))]
      temp=temp[,!grepl("conv",names(temp))]
      temp=temp[,!grepl("noboot",names(temp))]
      temp=temp[,!grepl("basic",names(temp))]
      
      res = rbind(res, temp)
      if (nrow(res) >= 10000) {break}
    }
    rm(i, temp)
    dim(res)
    
    
    options(scipen = 100, digits = 4)
    (dfl <- data.frame(
      Methode = c(    "Unadjusted",
                      "Elasticnet",
                      "Lasso",
                      "Neural network",
                      "Support vector machine",
                      "Super learner"
                      ,"Perfect"
                  
      ),
      
      MBp0 = c(calcul.mean.bias(parameter = res$mean.p0.unadjusted.bcv, estimand = mean(res_final$p0)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p0.elasticnet.bcv, estimand = mean(res_final$p0)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p0.lasso.bcv, estimand = mean(res_final$p0)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p0.mlp.bcv, estimand = mean(res_final$p0)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p0.svmRadial.bcv, estimand = mean(res_final$p0)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p0.sl.bcv[which(abs(res$mean.p0.sl.bcv) < 100)], estimand = mean(res_final$p0)) %>% round(2)
               ,calcul.mean.bias(parameter = res$mean.p0.perfect.bcv, estimand = mean(res_final$p0)) %>% round(2)
      ),
      
      MBp1 = c(calcul.mean.bias(parameter = res$mean.p1.unadjusted.bcv, estimand = mean(res_final$p1)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p1.elasticnet.bcv, estimand = mean(res_final$p1)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p1.lasso.bcv, estimand = mean(res_final$p1)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p1.mlp.bcv, estimand = mean(res_final$p1)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p1.svmRadial.bcv, estimand = mean(res_final$p1)) %>% round(2),
               calcul.mean.bias(parameter = res$mean.p1.sl.bcv[which(abs(res$mean.p1.sl.bcv) < 100)], estimand = mean(res_final$p1)) %>% round(2)
               ,calcul.mean.bias(parameter = res$mean.p1.perfect.bcv, estimand = mean(res_final$p1)) %>% round(2)
      ),
      
      logmOR = c((mean(res$mean.logOR.unadjusted.bcv) - mean(res_final$logOR)) %>% round(4),
                 (mean(res$mean.logOR.elasticnet.bcv) - mean(res_final$logOR)) %>% round(4),
                 (mean(res$mean.logOR.lasso.bcv) - mean(res_final$logOR)) %>% round(4),
                 (mean(res$mean.logOR.mlp.bcv) - mean(res_final$logOR)) %>% round(4),
                 (mean(res$mean.logOR.svmRadial.bcv) - mean(res_final$logOR)) %>% round(4),
                 (mean(res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)]) - mean(res_final$logOR)) %>% round(4)
                 ,(mean(res$mean.logOR.perfect.bcv) - mean(res_final$logOR)) %>% round(4)
      ),
      
      MBdelta = c(calcul.mean.bias(parameter = res$mean.delta.unadjusted.bcv, estimand = mean(res_final$delta)) %>% round(2),
                  calcul.mean.bias(parameter = res$mean.delta.elasticnet.bcv, estimand = mean(res_final$delta)) %>% round(2),
                  calcul.mean.bias(parameter = res$mean.delta.lasso.bcv, estimand = mean(res_final$delta)) %>% round(2),
                  calcul.mean.bias(parameter = res$mean.delta.mlp.bcv, estimand = mean(res_final$delta)) %>% round(2),
                  calcul.mean.bias(parameter = res$mean.delta.svmRadial.bcv, estimand = mean(res_final$delta)) %>% round(2),
                  calcul.mean.bias(parameter = res$mean.delta.sl.bcv[which(abs(res$mean.delta.sl.bcv) < 100)], estimand = mean(res_final$delta)) %>% round(2)
                  ,calcul.mean.bias(parameter = res$mean.delta.perfect.bcv, estimand = mean(res_final$delta)) %>% round(2)
      ),
      
      RMSE = c(calcul.root.mean.square.error(parameter = res$mean.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.root.mean.square.error(parameter = res$mean.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.root.mean.square.error(parameter = res$mean.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.root.mean.square.error(parameter = res$mean.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.root.mean.square.error(parameter = res$mean.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
               calcul.root.mean.square.error(parameter = res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
               ,calcul.root.mean.square.error(parameter = res$mean.logOR.perfect.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
      ),
      
      VEB = c(calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.unadjusted.bcv, sd.parameter = res$sd.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
              calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.elasticnet.bcv, sd.parameter = res$sd.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
              calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.lasso.bcv, sd.parameter = res$sd.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
              calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.mlp.bcv, sd.parameter = res$sd.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
              calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.svmRadial.bcv, sd.parameter = res$sd.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
              calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], sd.parameter = res$sd.logOR.sl.bcv[which(abs(res$mean.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
              ,calcul.variance.estimation.bias(mean.parameter = res$mean.logOR.perfect.bcv, sd.parameter = res$sd.logOR.perfect.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
      ),
      
      Coverage = c(calcul.coverage(borne.inf = res$inf.logOR.unadjusted.bcv, borne.sup = res$sup.logOR.unadjusted.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
                   calcul.coverage(borne.inf = res$inf.logOR.elasticnet.bcv, borne.sup = res$sup.logOR.elasticnet.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
                   calcul.coverage(borne.inf = res$inf.logOR.lasso.bcv, borne.sup = res$sup.logOR.lasso.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
                   calcul.coverage(borne.inf = res$inf.logOR.mlp.bcv, borne.sup = res$sup.logOR.mlp.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
                   calcul.coverage(borne.inf = res$inf.logOR.svmRadial.bcv, borne.sup = res$sup.logOR.svmRadial.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2),
                   calcul.coverage(borne.inf = res$inf.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], borne.sup = res$sup.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
                   ,calcul.coverage(borne.inf = res$inf.logOR.perfect.bcv, borne.sup = res$sup.logOR.perfect.bcv, estimand = mean(res_final$logOR[which(res_final$logOR != Inf)])) %>% round(2)
      ),
      
      Error = c(100 - calcul.puissance(borne.inf = res$inf.logOR.unadjusted.bcv, borne.sup = res$sup.logOR.unadjusted.bcv) %>% round(2),
                100 - calcul.puissance(borne.inf = res$inf.logOR.elasticnet.bcv, borne.sup = res$sup.logOR.elasticnet.bcv) %>% round(2),
                100 - calcul.puissance(borne.inf = res$inf.logOR.lasso.bcv, borne.sup = res$sup.logOR.lasso.bcv) %>% round(2),
                100 - calcul.puissance(borne.inf = res$inf.logOR.mlp.bcv, borne.sup = res$sup.logOR.mlp.bcv) %>% round(2),
                100 - calcul.puissance(borne.inf = res$inf.logOR.svmRadial.bcv, borne.sup = res$sup.logOR.svmRadial.bcv) %>% round(2),
                100 - calcul.puissance(borne.inf = res$inf.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)], borne.sup = res$sup.logOR.sl.bcv[which(abs(res$inf.logOR.sl.bcv) < 100)]) %>% round(2)
                ,100 - calcul.puissance(borne.inf = res$inf.logOR.perfect.bcv, borne.sup = res$sup.logOR.perfect.bcv) %>% round(2)
      ),
      
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
              ,calcul.rss(mean.parameter = res$mean.logOR.perfect.bcv, sd.parameter = res$sd.logOR.perfect.bcv,
                          mean.parameter.unadjusted = res$mean.logOR.unadjusted.bcv, sd.parameter.unadjusted = res$sd.logOR.unadjusted.bcv) %>% round(2)
      )
      
      
    ))
    dfl
    dfl[-nrow(dfl),] #Without the last line of the perfect model (not in paper)
    
    
    
    
    dfl[,2]=paste0(sprintf("%.2f", dfl[,2]), "%")
    dfl[,3]=paste0(sprintf("%.2f", dfl[,3]), "%")
    dfl[,4]=sprintf("%.4f", dfl[,4])
    dfl[,5]=paste0(sprintf("%.2f", dfl[,5]), "%")
    dfl[,7]=paste0(sprintf("%.2f", dfl[,7]), "%")
    dfl[,8]=paste0(sprintf("%.2f", dfl[,8]), "%")
    dfl[,9]=paste0(sprintf("%.2f", dfl[,9]), "%")
    dfl[,10]=paste0(sprintf("%.2f", dfl[,10]), "%")
    
    
    dfl[,1]=c("Unadjusted","Elasticnet","Lasso","Neural network","Support vector machine","Super learner","Perfect")
    
    print(xtable(dfl[-nrow(dfl),], type = "latex"), file = paste0(paste0("Latex_complex_n", N.effectif, "_ate", Size.Effect,".tex")),include.rownames=FALSE)
    
    
    
  }
}



