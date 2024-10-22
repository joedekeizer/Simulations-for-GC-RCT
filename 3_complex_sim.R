################################################################################
###### (3) Applying the different methods on the simulated data
################################################################################
### This code will create a folder in the specified work directory named /results_complex_nXX_ateYY where XX is the sample size (N.effectif) and YY the effect size (Size.Effect)
### In the folder it will create a .csv file with the results of all the employed methods on the simulated data
### Paralellization is highly recommended to reduce the processing time, the code is written so that it saves the progress as it runs with the possibility to pause and resume (skip data already with results)
################################################################################
library(dplyr)
library(MASS)
library(splines)
library(caret)
library(SuperLearner)
library(glmnet)
library(doParallel)
library(cvAUC)

### Initialization of the parameters
SizeEffect <- "log(3)" #The size effect (character)
N.effectif <- 200 #The sample size (numeric)
N.stop <- 10000 #The number of simulated data sets (numeric)
B <- 500 #The number of bootstrap samples
V.crossvalid <- 20 #The number of cross-validations
path <- "/home/Simulations/" #Work directory


### Setting the work directory to the folder with the simulated datasets
Size.Effect <- gsub("[[:punct:]]", "", SizeEffect) #Removing any special character for folder name
setwd(paste0(path, "/complex_n", N.effectif, "_ate", Size.Effect, "/")) #set to work directory to the location of the simulated data


### Create the folder for the results if it does not already exist in the work directory
pathresults <- paste0(path, "results_complex_n", N.effectif, "_ate", Size.Effect, "/")
ifelse(!dir.exists(file.path(pathresults)), dir.create(file.path(pathresults)), FALSE)


### Parallelisation of the code (PSOCK if Windows and FORK if Linux environment)
nb.cluster <- parallel::detectCores() - 1 #All cores minus one
if(.Platform[[1]] == "windows") {cl <- makeCluster(nb.cluster, type="PSOCK") } else{ cl <- makeCluster(nb.cluster, type="FORK") }
registerDoParallel(cl)
if(.Platform[[1]] == "windows") {clusterEvalQ(cl, {library(dplyr);library(MASS);library(splines);library(caret);library(SuperLearner);library(glmnet);library(cvAUC)})}


################################################################################

### Function to use our ML methods on the datasets and save the results of all 500 bootstraps 
iteration <- function(iter) 
{
  
  ### Functions used for reporting the results and for the Superlearner 
  calcul.logOR <- function(p0, p1)
  {
    return(log((p1*(1-p0))/(p0*(1-p1))))
  }
  
  ### SL elasticnet
  SL.elasticnet.caret <- function (Y, X, newX, family,  ...)
  {
    X <- model.matrix(~ -1 + ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                      x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                      bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                      bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)), X)
    
    newX <- model.matrix(~ -1 + ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                         x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                         bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                         bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)),  newX)
    
    fitElastic <- glmnet(x = X, y = Y, family = "binomial", alpha = elasticnet.param$bestTune$alpha, lambda = elasticnet.param$bestTune$lambda, penalty.factor = c(0, rep(1, 78)) )
    
    pred <- predict(fitElastic, newx = newX, type = "response")
    fit <- list(object = fitElastic)
    class(fit) <- "SL.elasticnet.caret"
    out <- list(pred = pred, fit = fit)
    return(out)
  }
  
  ### SL lasso
  SL.lasso.caret <- function (Y, X, newX, family,  ...)
  {
    X <- model.matrix(~ -1 + ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                      x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                      bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                      bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)), X)
    
    newX <- model.matrix(~ -1 + ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                         x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                         bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                         bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)),  newX)
    
    fitLasso <- glmnet(x = X, y = Y, family = "binomial", alpha = 1,  lambda = lasso.param$bestTune$lambda,  penalty.factor = c(0, rep(1, 78)))
    
    pred <- predict(fitLasso, newx = newX, type = "response")
    fit <- list(object = fitLasso)
    class(fit) <- "SL.lasso.caret"
    out <- list(pred = pred, fit = fit)
    return(out)
  }
  
  ### SL neural network
  SL.mlp.caret <- function (Y, X, newX, family, ...)
  {
    fit.mlp <- RSNNS::mlp(x = X,  y = Y, size = mlp.param$bestTune$size, trace = FALSE)
    
    pred <- predict(fit.mlp, newdata = newX, type = "raw")
    fit <- list(object = fit.mlp)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.mlp.caret")
    return(out)
  }
  
  ### SL support vector machine
  SL.ksvm.caret <- function (Y, X, newX, family, type = "C-svc", kernel = "rbfdot", kpar = "automatic", nu = 0.2, epsilon = 0.1, cross = 0, prob.model = family$family == "binomial", class.weights = NULL, cache = 40, tol = 0.001, shrinking = T, ...)
  {
    if (!is.matrix(X)) {X <- model.matrix( ~ ., data = X); X <- X[,-1]}
    Y <- as.factor(Y)
    predict_type <- "probabilities"
    model <- kernlab::ksvm(X, Y, type = "C-svc", kernel = "rbfdot",  kpar = list(sigma = svmRadial.param$bestTune$sigma),
                           C = svmRadial.param$bestTune$C, nu = nu, epsilon = epsilon, prob.model = prob.model, class.weights = class.weights)
    
    if (!is.matrix(newX)) {newX <-  model.matrix( ~ ., data = newX); newX <- newX[,-1, drop = FALSE]}
    
    pred <- kernlab::predict(model, newX, predict_type)[, 2]
    fit <- list(object = model, family = family)
    out <- list(pred = pred, fit = fit)
    class(out$fit) = "SL.ksvm.caret"
    return(out)
  }
  
  ### Load a simulated dataset
  load(liste[iter])
  data.obs1 <- data.obs0 <- data.obs
  
  ### Creating the counterfactual datasets for the GC
  data.obs1$ttt <- 1
  data.obs0$ttt <- 0
  
  ### Preparation for the ML
  base.train <- data.obs
  base.train$outcome[base.train$outcome == 1] <- "yes"
  base.train$outcome[base.train$outcome == 0] <- "no"
  base.train$outcome <- as.factor(base.train$outcome)
  
  ### Control parameters for the predictions
  control <- trainControl(allowParallel = FALSE,
                          verboseIter = FALSE,
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          method = "cv")
  
  ### Tuning parameters by maximizing the AUC of the ROC curve with 20-fold cross-validation
  mlp.param <- train(outcome ~ ., data = base.train, method = 'mlp', tuneGrid = expand.grid(size=1:20), metric = "ROC", trControl = control)
  elasticnet.param <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                              x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                              bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                              bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)), data = base.train,
                           method = 'glmnet',  tuneLength = 20, metric = "ROC", trControl = control,
                           family = "binomial", penalty.factor = c(0, rep(1, 78)))
  lasso.param <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                          x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                          bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                          bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)),
                       data = base.train, method = 'glmnet',   tuneGrid = expand.grid(.alpha = 1, .lambda = unique(elasticnet.param$results$lambda)),
                       metric = "ROC", trControl = control, family = "binomial", penalty.factor = c(0, rep(1, 78)) )
  svmRadial.param <- train(outcome ~ .,  data = base.train, method = 'svmRadialSigma', trace = FALSE, tuneLength = 20,
                            metric = "ROC", trControl = control,  allowParallel = FALSE)

  y.sl <- as.numeric(base.train$outcome) - 1
  x.sl <- base.train[, names(base.train) != "outcome"]
  x.sl.new <- data.frame(rbind(x.sl, x.sl))
  x.sl.new$ttt[1:N.effectif] <- 1 #1
  x.sl.new$ttt[(N.effectif + 1):(2 * N.effectif)] <- 0 #0

  sl <- SuperLearner(Y = y.sl, X = x.sl,  newX = x.sl.new, family = binomial(),
                     SL.library = list("SL.mlp.caret", "SL.elasticnet.caret",
                                       "SL.lasso.caret", "SL.ksvm.caret"),
                     method = method.AUC(optim_method = "BFGS", bounds = c(-Inf, Inf)),
                     cvControl = list(V = V.crossvalid))
  

  ### Initialization of the results
  p0.unadjusted.bcv <- p1.unadjusted.bcv <- logOR.unadjusted.bcv <- delta.unadjusted.bcv <-
  p0.elasticnet.bcv <- p1.elasticnet.bcv <- logOR.elasticnet.bcv <- delta.elasticnet.bcv <-
  p0.lasso.bcv <- p1.lasso.bcv <- logOR.lasso.bcv <- delta.lasso.bcv <-
  p0.mlp.bcv <- p1.mlp.bcv <- logOR.mlp.bcv <- delta.mlp.bcv <-
  p0.svmRadial.bcv <- p1.svmRadial.bcv <- logOR.svmRadial.bcv <- delta.svmRadial.bcv <-
  p0.sl.bcv <- p1.sl.bcv <- logOR.sl.bcv <- delta.sl.bcv <-
  NA
  
  
  ### Bootstrapping on the dataset
  for(b in 1:B)
  {
    ### Create the models on the learn sample
    id <- sample(1:N.effectif, size = N.effectif, replace = TRUE)
    learn <- data.obs[id,]
    
    ### Predict the parameters on the remaining observations
    valid <- data.obs[-sort(unique(id)),]
    valid0 <- valid1 <- valid
    valid0$ttt <- 0
    valid1$ttt <- 1
    
    ### Estimate unadjusted p0 et p1
    mod.unadjusted <- glm(outcome ~ ttt, data = learn, family = binomial(link = "logit"))
    
    p0.unadjusted.bcv[b] <- predict(mod.unadjusted, newdata = valid0, type = "response") %>% mean
    p1.unadjusted.bcv[b] <- predict(mod.unadjusted, newdata = valid1, type = "response") %>% mean
    logOR.unadjusted.bcv[b] <- calcul.logOR(p0.unadjusted.bcv[b], p1.unadjusted.bcv[b])
    delta.unadjusted.bcv[b] <- p1.unadjusted.bcv[b] - p0.unadjusted.bcv[b]
    
    ### Preparation for the ML
    learn.caret <- learn
    learn.caret$outcome[learn.caret$outcome == 1] <- "yes"
    learn.caret$outcome[learn.caret$outcome == 0] <- "no"
    learn.caret$outcome <- as.factor(learn.caret$outcome)

    valid.caret <- valid
    valid.caret$outcome[valid.caret$outcome == 1] <- "yes"
    valid.caret$outcome[valid.caret$outcome == 0] <- "no"
    valid.caret$outcome <- as.factor(valid.caret$outcome)
    valid.caret0 <- valid.caret1 <- valid.caret
    valid.caret0$ttt <- 0
    valid.caret1$ttt <- 1
    
    control <- trainControl(allowParallel = FALSE,
                            verboseIter = FALSE,
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary,
                            method = "cv")
    
    
    ### Estimate via elasticnet (using interactions and splines)
    mod.elasticnet <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                               x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                               bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                               bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)), data = learn.caret,
                            method = 'glmnet',
                            tuneGrid = expand.grid(alpha = elasticnet.param$bestTune$alpha, lambda = elasticnet.param$bestTune$lambda),
                            metric = "ROC", trControl = control,
                            family = "binomial", penalty.factor = c(0, rep(1, 78)))

    p0.elasticnet.bcv[b] <- mean(predict(object = mod.elasticnet, newdata = valid.caret0,  type = "prob")[,2])
    p1.elasticnet.bcv[b] <- mean(predict(object = mod.elasticnet, newdata = valid.caret1,  type = "prob")[,2])
    logOR.elasticnet.bcv[b] <- calcul.logOR(p0.elasticnet.bcv[b], p1.elasticnet.bcv[b])
    delta.elasticnet.bcv[b] <- p1.elasticnet.bcv[b] - p0.elasticnet.bcv[b]
    
    
    ### Estimate via mlp
    mod.mlp <-  train(outcome ~ ., data = learn.caret,
                      method = 'mlp',
                      tuneGrid = expand.grid(size = mlp.param$bestTune$size),
                      metric = "ROC", trControl = control)

    p0.mlp.bcv[b] <- mean(predict(object = mod.mlp, newdata = valid.caret0,  type = "prob")[,2])
    p1.mlp.bcv[b] <- mean(predict(object = mod.mlp, newdata = valid.caret1,  type = "prob")[,2])
    logOR.mlp.bcv[b] <- calcul.logOR(p0.mlp.bcv[b], p1.mlp.bcv[b])
    delta.mlp.bcv[b] <- p1.mlp.bcv[b] - p0.mlp.bcv[b]
    
    
    ### Estimate via Lasso (using interactions and splines)
    mod.lasso <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
                                          x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
                                          bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
                                          bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)),
                       data = learn.caret, method = 'glmnet',   tuneGrid = expand.grid(alpha=lasso.param$bestTune$alpha, lambda=lasso.param$bestTune$lambda),
                       metric = "ROC", trControl = control, family = "binomial", penalty.factor = c(0, rep(1, 78)) )

    p0.lasso.bcv[b] <- mean(predict(object = mod.lasso, newdata = valid.caret0,  type = "prob")[,2])
    p1.lasso.bcv[b] <- mean(predict(object = mod.lasso, newdata = valid.caret1,  type = "prob")[,2])
    logOR.lasso.bcv[b] <- calcul.logOR(p0.lasso.bcv[b], p1.lasso.bcv[b])
    delta.lasso.bcv[b] <- p1.lasso.bcv[b] - p0.lasso.bcv[b]
    
    
    ### Estimate via svm
    mod.svmRadial <-  train(outcome ~ .,  data = learn.caret, method = 'svmRadialSigma', trace = FALSE,
                            tuneGrid = expand.grid(sigma=svmRadial.param$bestTune$sigma, C = svmRadial.param$bestTune$C),
                            metric = "ROC", trControl = control,  allowParallel = FALSE)

    p0.svmRadial.bcv[b] <- mean(predict(object = mod.svmRadial, newdata = valid.caret0,  type = "prob")[,2])
    p1.svmRadial.bcv[b] <- mean(predict(object = mod.svmRadial, newdata = valid.caret1,  type = "prob")[,2])
    logOR.svmRadial.bcv[b] <- calcul.logOR(p0.svmRadial.bcv[b], p1.svmRadial.bcv[b])
    delta.svmRadial.bcv[b] <- p1.svmRadial.bcv[b] - p0.svmRadial.bcv[b]
    

    ### Estimate via Superlearner
    outcome.sl <-  as.numeric(learn.caret$outcome) - 1
    covariables.sl <-  learn.caret[, names(learn.caret) != "outcome"]

    N.valid <- length(valid.caret$outcome)
    outcome.sl.valid <-  as.numeric(valid.caret$outcome) - 1
    covariables.sl.valid <-  valid.caret[, names(valid.caret) != "outcome"]
    covariables.sl.valid <- data.frame(rbind(covariables.sl.valid, covariables.sl.valid))
    covariables.sl.valid$ttt[1:N.valid] <- 1
    covariables.sl.valid$ttt[(N.valid + 1):(2 * N.valid)] <- 0
    
    sl <- SuperLearner(Y = outcome.sl, X = covariables.sl,  newX = covariables.sl.valid,
                       family = binomial(),
                       SL.library = list("SL.mlp.caret", "SL.elasticnet.caret", "SL.lasso.caret", "SL.ksvm.caret"),
                       method = method.AUC(optim_method = "BFGS", bounds = c(-Inf, Inf)),
                       cvControl = list(V = V.crossvalid))
    pred.sl.bcv <- predict(sl, onlySL = TRUE)$pred
    p0.sl.bcv[b] <- mean(pred.sl.bcv[(N.valid + 1):(2 * N.valid)])
    p1.sl.bcv[b] <- mean(pred.sl.bcv[1:N.valid])
    logOR.sl.bcv[b] <- calcul.logOR(p0 = p0.sl.bcv[b], p1 = p1.sl.bcv[b])
    delta.sl.bcv[b] <- p1.sl.bcv[b] - p0.sl.bcv[b]
    
  }
  
  ### Saving the results of each bootstrap
  res.iter <- data.frame(
    iteration = iter,
    n0 = sum(1*I(data.obs$ttt == 0)),
    n1 = sum(1*I(data.obs$ttt == 1)),
    
    b.unadjusted.bcv = sum(1*I(!is.na(logOR.unadjusted.bcv))),
    mean.p0.unadjusted.bcv = mean(p0.unadjusted.bcv, na.rm = T),
    mean.p1.unadjusted.bcv = mean(p1.unadjusted.bcv, na.rm = T),
    mean.logOR.unadjusted.bcv = mean(logOR.unadjusted.bcv, na.rm = T),
    mean.delta.unadjusted.bcv = mean(delta.unadjusted.bcv, na.rm = T),
    sd.p0.unadjusted.bcv = sd(p0.unadjusted.bcv, na.rm = T),
    sd.p1.unadjusted.bcv = sd(p1.unadjusted.bcv, na.rm = T),
    sd.logOR.unadjusted.bcv = sd(logOR.unadjusted.bcv, na.rm = T),
    sd.delta.unadjusted.bcv = sd(delta.unadjusted.bcv, na.rm = T),
    inf.p0.unadjusted.bcv = quantile(p0.unadjusted.bcv, probs = 0.025, na.rm = T),
    inf.p1.unadjusted.bcv = quantile(p1.unadjusted.bcv, probs = 0.025, na.rm = T),
    inf.logOR.unadjusted.bcv = quantile(logOR.unadjusted.bcv, probs = 0.025, na.rm = T),
    inf.delta.unadjusted.bcv = quantile(delta.unadjusted.bcv, probs = 0.025, na.rm = T),
    sup.p0.unadjusted.bcv = quantile(p0.unadjusted.bcv, probs = 0.975, na.rm = T),
    sup.p1.unadjusted.bcv = quantile(p1.unadjusted.bcv, probs = 0.975, na.rm = T),
    sup.logOR.unadjusted.bcv = quantile(logOR.unadjusted.bcv, probs = 0.975, na.rm = T),
    sup.delta.unadjusted.bcv = quantile(delta.unadjusted.bcv, probs = 0.975, na.rm = T),
    
    b.elasticnet.bcv = sum(1*I(!is.na(logOR.elasticnet.bcv))),
    mean.p0.elasticnet.bcv = mean(p0.elasticnet.bcv, na.rm = T),
    mean.p1.elasticnet.bcv = mean(p1.elasticnet.bcv, na.rm = T),
    mean.logOR.elasticnet.bcv = mean(logOR.elasticnet.bcv, na.rm = T),
    mean.delta.elasticnet.bcv = mean(delta.elasticnet.bcv, na.rm = T),
    sd.p0.elasticnet.bcv = sd(p0.elasticnet.bcv, na.rm = T),
    sd.p1.elasticnet.bcv = sd(p1.elasticnet.bcv, na.rm = T),
    sd.logOR.elasticnet.bcv = sd(logOR.elasticnet.bcv, na.rm = T),
    sd.delta.elasticnet.bcv = sd(delta.elasticnet.bcv, na.rm = T),
    inf.p0.elasticnet.bcv = quantile(p0.elasticnet.bcv, probs = 0.025, na.rm = T),
    inf.p1.elasticnet.bcv = quantile(p1.elasticnet.bcv, probs = 0.025, na.rm = T),
    inf.logOR.elasticnet.bcv = quantile(logOR.elasticnet.bcv, probs = 0.025, na.rm = T),
    inf.delta.elasticnet.bcv = quantile(delta.elasticnet.bcv, probs = 0.025, na.rm = T),
    sup.p0.elasticnet.bcv = quantile(p0.elasticnet.bcv, probs = 0.975, na.rm = T),
    sup.p1.elasticnet.bcv = quantile(p1.elasticnet.bcv, probs = 0.975, na.rm = T),
    sup.logOR.elasticnet.bcv = quantile(logOR.elasticnet.bcv, probs = 0.975, na.rm = T),
    sup.delta.elasticnet.bcv = quantile(delta.elasticnet.bcv, probs = 0.975, na.rm = T),
    
    b.mlp.bcv = sum(1*I(!is.na(logOR.mlp.bcv))),
    mean.p0.mlp.bcv = mean(p0.mlp.bcv, na.rm = T),
    mean.p1.mlp.bcv = mean(p1.mlp.bcv, na.rm = T),
    mean.logOR.mlp.bcv = mean(logOR.mlp.bcv, na.rm = T),
    mean.delta.mlp.bcv = mean(delta.mlp.bcv, na.rm = T),
    sd.p0.mlp.bcv = sd(p0.mlp.bcv, na.rm = T),
    sd.p1.mlp.bcv = sd(p1.mlp.bcv, na.rm = T),
    sd.logOR.mlp.bcv = sd(logOR.mlp.bcv, na.rm = T),
    sd.delta.mlp.bcv = sd(delta.mlp.bcv, na.rm = T),
    inf.p0.mlp.bcv = quantile(p0.mlp.bcv, probs = 0.025, na.rm = T),
    inf.p1.mlp.bcv = quantile(p1.mlp.bcv, probs = 0.025, na.rm = T),
    inf.logOR.mlp.bcv = quantile(logOR.mlp.bcv, probs = 0.025, na.rm = T),
    inf.delta.mlp.bcv = quantile(delta.mlp.bcv, probs = 0.025, na.rm = T),
    sup.p0.mlp.bcv = quantile(p0.mlp.bcv, probs = 0.975, na.rm = T),
    sup.p1.mlp.bcv = quantile(p1.mlp.bcv, probs = 0.975, na.rm = T),
    sup.logOR.mlp.bcv = quantile(logOR.mlp.bcv, probs = 0.975, na.rm = T),
    sup.delta.mlp.bcv = quantile(delta.mlp.bcv, probs = 0.975, na.rm = T),
    
    b.lasso.bcv = sum(1*I(!is.na(logOR.lasso.bcv))),
    mean.p0.lasso.bcv = mean(p0.lasso.bcv, na.rm = T),
    mean.p1.lasso.bcv = mean(p1.lasso.bcv, na.rm = T),
    mean.logOR.lasso.bcv = mean(logOR.lasso.bcv, na.rm = T),
    mean.delta.lasso.bcv = mean(delta.lasso.bcv, na.rm = T),
    sd.p0.lasso.bcv = sd(p0.lasso.bcv, na.rm = T),
    sd.p1.lasso.bcv = sd(p1.lasso.bcv, na.rm = T),
    sd.logOR.lasso.bcv = sd(logOR.lasso.bcv, na.rm = T),
    sd.delta.lasso.bcv = sd(delta.lasso.bcv, na.rm = T),
    inf.p0.lasso.bcv = quantile(p0.lasso.bcv, probs = 0.025, na.rm = T),
    inf.p1.lasso.bcv = quantile(p1.lasso.bcv, probs = 0.025, na.rm = T),
    inf.logOR.lasso.bcv = quantile(logOR.lasso.bcv, probs = 0.025, na.rm = T),
    inf.delta.lasso.bcv = quantile(delta.lasso.bcv, probs = 0.025, na.rm = T),
    sup.p0.lasso.bcv = quantile(p0.lasso.bcv, probs = 0.975, na.rm = T),
    sup.p1.lasso.bcv = quantile(p1.lasso.bcv, probs = 0.975, na.rm = T),
    sup.logOR.lasso.bcv = quantile(logOR.lasso.bcv, probs = 0.975, na.rm = T),
    sup.delta.lasso.bcv = quantile(delta.lasso.bcv, probs = 0.975, na.rm = T),
    
    b.svmRadial.bcv = sum(1*I(!is.na(logOR.svmRadial.bcv))),
    mean.p0.svmRadial.bcv = mean(p0.svmRadial.bcv, na.rm = T),
    mean.p1.svmRadial.bcv = mean(p1.svmRadial.bcv, na.rm = T),
    mean.logOR.svmRadial.bcv = mean(logOR.svmRadial.bcv, na.rm = T),
    mean.delta.svmRadial.bcv = mean(delta.svmRadial.bcv, na.rm = T),
    sd.p0.svmRadial.bcv = sd(p0.svmRadial.bcv, na.rm = T),
    sd.p1.svmRadial.bcv = sd(p1.svmRadial.bcv, na.rm = T),
    sd.logOR.svmRadial.bcv = sd(logOR.svmRadial.bcv, na.rm = T),
    sd.delta.svmRadial.bcv = sd(delta.svmRadial.bcv, na.rm = T),
    inf.p0.svmRadial.bcv = quantile(p0.svmRadial.bcv, probs = 0.025, na.rm = T),
    inf.p1.svmRadial.bcv = quantile(p1.svmRadial.bcv, probs = 0.025, na.rm = T),
    inf.logOR.svmRadial.bcv = quantile(logOR.svmRadial.bcv, probs = 0.025, na.rm = T),
    inf.delta.svmRadial.bcv = quantile(delta.svmRadial.bcv, probs = 0.025, na.rm = T),
    sup.p0.svmRadial.bcv = quantile(p0.svmRadial.bcv, probs = 0.975, na.rm = T),
    sup.p1.svmRadial.bcv = quantile(p1.svmRadial.bcv, probs = 0.975, na.rm = T),
    sup.logOR.svmRadial.bcv = quantile(logOR.svmRadial.bcv, probs = 0.975, na.rm = T),
    sup.delta.svmRadial.bcv = quantile(delta.svmRadial.bcv, probs = 0.975, na.rm = T),

    b.sl.bcv = sum(1*I(!is.na(logOR.sl.bcv))),
    mean.p0.sl.bcv = mean(p0.sl.bcv, na.rm = T),
    mean.p1.sl.bcv = mean(p1.sl.bcv, na.rm = T),
    mean.logOR.sl.bcv = mean(logOR.sl.bcv, na.rm = T),
    mean.delta.sl.bcv = mean(delta.sl.bcv, na.rm = T),
    sd.p0.sl.bcv = sd(p0.sl.bcv, na.rm = T),
    sd.p1.sl.bcv = sd(p1.sl.bcv, na.rm = T),
    sd.logOR.sl.bcv = sd(logOR.sl.bcv, na.rm = T),
    sd.delta.sl.bcv = sd(delta.sl.bcv, na.rm = T),
    inf.p0.sl.bcv = quantile(p0.sl.bcv, probs = 0.025, na.rm = T),
    inf.p1.sl.bcv = quantile(p1.sl.bcv, probs = 0.025, na.rm = T),
    inf.logOR.sl.bcv = quantile(logOR.sl.bcv, probs = 0.025, na.rm = T),
    inf.delta.sl.bcv = quantile(delta.sl.bcv, probs = 0.025, na.rm = T),
    sup.p0.sl.bcv = quantile(p0.sl.bcv, probs = 0.975, na.rm = T),
    sup.p1.sl.bcv = quantile(p1.sl.bcv, probs = 0.975, na.rm = T),
    sup.logOR.sl.bcv = quantile(logOR.sl.bcv, probs = 0.975, na.rm = T),
    sup.delta.sl.bcv = quantile(delta.sl.bcv, probs = 0.975, na.rm = T)
    
  )
  
  ### Saving the results of all bootstraps for one dataset into a .csv file
  write.table(res.iter, paste0(pathresults, liste[iter], ".csv"), sep = ";", row.names = F)
  return(res.iter)
}


################################################################################

### List all the simulated datasets
liste <- list.files() %>% as.character #list all 10000 simulated datasets

### Skipping datasets where the results have already been created
totlist <- 1:N.stop
listmissing <- totlist[!totlist %in% as.numeric(gsub(".Rdata.csv","",gsub(paste0("sim_n",N.effectif,"_ate",Size.Effect,"_"),"",list.files(pathresults))))]


### Running the function with parallel processing (highly recommended)
if(.Platform[[1]] == "windows") {
  ### Windows
  foreach(i = listmissing, .combine = rbind, .inorder = FALSE, .verbose = T,
          .export=ls(envir=globalenv()),
          .packages = c("dplyr", "MASS", "splines", "caret", "SuperLearner", "glmnet", "cvAUC")
  ) %dopar% {iteration(i)}
} else {
  ### Linux
  foreach(i = listmissing, .combine = rbind, .inorder = FALSE, .verbose = T
  ) %dopar% {iteration(i)}
}

stopCluster(cl)
