R codes used in the article titled: G-computation for increasing performances of randomized clinical trials with binary response
================


## Description

In the article G-computation for increasing performances of clinical trials with individual randomization and binary response (DOI. XX), we investigated by simulations the performances of G-computation according to various machine learning algorithms in the context of randomized clinical trials with a 1:1 allocation ratio, a superiority design and binary outcomes. 
This repository contains all the R code used in the article, organized into four different steps and R scripts. The code is optimized for the Linux environment but is also executable on Windows. To enhance processing speed, parallelization using the `doParallel` package is proposed.
The provided R scripts contain the code for the complex scenario with n = 200 and mOR = 1.9. At the beginning of each R script, a list of parameters is available for simulating the other complex scenarios shown in the article.
All analyses were performed using R version 4.3.0.


## Data generation

Two distinct scenarios were simulated. In the provided code, we consider the complex scenario where the estimation of the outcome model was not straightforward. In the article, we considered different marginal treatment effects and sample sizes. For all of the scenarios 10,000 data sets were generated.

![dagcomplexe](https://github.com/user-attachments/assets/82e7b683-0814-45fd-8ca2-1b675d03dfc0)

The first R script, named `1_sim_data.R`, contains the code to simulate data sets with the relationships depicted in the DAG of the complex scenario. In this script, only the first few lines need to be modified to run the code with different effect and sample sizes.
When executed, this script will create a folder named `/complex_nXX_ateYY`, where `XX` is the sample size (`N.effectif`) and `YY` the effect size (`Size.Effect`), containing all the simulated data sets in `.Rdata` format.
To simulate the data sets of the simple scenario, only the `sim.base` function needs to be modified with the variables and coefficients described in Supplementary table A2 of the article (see tables at the bottom of the page).

``` r
### Initialization
path <- "/home/Simulations/" #Work directory
Size.Effect <- "log3" #The beta4 coefficient from Supplementary Table A1 for folder names (character format)
sizeEffect <- log(3) #Same beta4 coefficient but numeric for the simulations (numeric)
N.effectif <- 200 #Sample size (numeric)
N.stop <- 10000 #10k data sets to simulate and save as Rdata (data.obs)
```


## Estimating the theoretical estimates of the Data

In the file named `2_estimate_param_ATE.R`, we provide a code to obtain the theoretical values of the marginal effects by the mean of the unadjusted estimations for 1,000,000 simulated data sets. As in the previous script, only the parameters at the beginning need to be specified to obtain the outcomes for different scenarios.
Similar to the earlier script, the code contains the causal relationships between the treatment arm, the covariates and outcome for the complex scenario. To obtain the estimates for the simple scenario, only the `sim.base` function needs to be modified.
This script will generate an `.Rdata` file with 1,000,000 rows, each containing the number of simulated observations in each arm, the proportions of the simulated outcomes per arm, the log OR and the difference in proportions.


## Results for the different Machine learning techniques

We considered several models and algorithms (learners) to fit the outcome model. All analyses were performed using R version 4.3.0, using the `caret` package with a tuning grid of length equal to 20. Below an overview of the learners used:

• Lasso logistic regression. L1 regularization allows for the selection of the predictors. To establish a flexible model, we considered all possible interactions between the treatment arm A and the covariates X. Additionally, we used B-splines for the continuous covariates. The `glmnet` package was used. The penalization of the L1 norm was the only tuning parameter.
• Elasticnet logistic regression. This approach mirrors the logistic regression mentioned earlier but incorporates both the L1 and the L2 regularizations.
• Neural network. We chose one hidden layer, which represents one of the most common network architectures. Its size constitutes the tuning parameter. The `nnet` package was used.
• Support vector machine. To relax the linear assumption, we opted for the radial basis function kernel. The `svmRadial` function of the `kernlab` package was used. It requires two tuning parameters: the cost penalty of miss-classification and the flexibility of the classification.
• Super learner. We also tested a super learner with the ensemble of the previous ML techniques. Super learner consisted in a weighting average of the learner-specific predictions by using a weighted linear predictor. In alignment with our previous choices, we estimated the weights by maximizing the average AUC through a 20-fold cross-validation. We used the `SuperLearner` package.

The third script, `3_complex_sim.R`, applies all the previously described techniques to the 10,000 simulated data sets generated using the first script (1_sim_data.R). The script then creates a new folder named `/results_complex_nXX_ateYY` where `XX` is the sample size (`N.effectif`) and `YY` the effect size (`Size.Effect`), with all the results of the techniques applied to each simulated data set. The output is a `.csv file` containing a list of estimates for each technique and for for each simulated data set (as shown below).

``` r
### Results of the Elasticnet technique (500 bootstrap cross-validation)
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
sup.delta.elasticnet.bcv = quantile(delta.elasticnet.bcv, probs = 0.975, na.rm = T)
```

Note that, for both penalized methods, we considered all possible interactions between the treatment arm and the covariates. Additionally, we used B-splines for the continuous covariates, with three degrees of freedom.
To use this script on the simple scenario, the model specifications for both penalized methods need to be adjusted accordingly. This modification is needed in the function for the Superlearner, for the tuning parameters and inside the bootstrapping process to predict the estimates. See R code below with an example for obtaining the tuning parameters for the Elasticnet model in both the complex and the simple scenario. Note not only the difference in the formula but also the `penality.factor`, which corresponds to the total number of covariates. Un this case, the splines had three degrees of freedom.

``` r
### Tuning parameters by maximizing the AUC of the ROC curve with 20-fold cross-validation
### Elasticnet for the complex scenario :
elasticnet.param <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + bs(x5, df = 3) +
						x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) +
						bs(x11, df = 3) + x12 + bs(x14, df = 3) + bs(x15, df = 3) +
						bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)), data = base.train,
				method = 'glmnet',  tuneLength = 20, metric = "ROC", trControl = control,
				family = "binomial", penalty.factor = c(0, rep(1, 78)))

### Elasticnet for the simple scenario :
elasticnet.param = train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) + x4 + x5 + x6), data = base.train,
				method = 'glmnet',  tuneLength = 20, metric = "ROC", trControl = control,
				family = "binomial", penalty.factor = c(0,rep(1,24)))
``` r			   

This R script has a long processing time, so parallel processing is highly recommended. The code was written to allow for stopping and resuming execution from where it was left off.


## Summarizing all the Results
The fourth R script, named `4_read_sims.R`, contains all the functions necessary to summarize the different performance criteria used for each of the techniques on the 10,000 simulated data sets. It creates a table as a data frame in R, with the unadjusted estimates and the different techniques represented in the rows, and the different criteria in the columns.
Finally, tt then saves the table as a `.tex` file.


## Models used for the simulations in the complex and the simple scenario
![tablea1](https://github.com/user-attachments/assets/2edbb575-5454-4c73-a0d9-ed14c1ba4e6b)

![tablea2](https://github.com/user-attachments/assets/b02bdc24-c133-4185-affb-a8b112e12869)