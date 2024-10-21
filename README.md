R scripts for data simulation in de Keizer et al. (submitted papier, 2024)
================


## Description

We aimed to investigate by simulations the performances of G-computation according to various machine learning algorithms in the context of randomized clinical trials with a 1:1 allocation ratio, a superiority design and binary outcomes. 


## Simulating the 10,000 datasets

This repository contains all the R scripts for the complex scenario with sample size $n$ = 200 and marginal $mOR$ = 1.9. The code is optimized for the Linux environment but is also executable on Windows.  At the beginning of each R script, a list of parameters is available for simulating the other complex scenarios shown in the article. The first file `1_sim_data.R` contains the scripts to simulate data sets with the following complex scenario. 

<p align="center">
<img src="./pictures/DAGcomplexe.png" alt="image" width="700" height="auto">
</p>

Only the first lines need to be modified to run the code with different effect and sample sizes. When executed, this script will create a folder named `/complex_nXX_ateYY`, where `XX` is the sample size (`N.effectif`) and `YY` the effect size (`Size.Effect`), containing all the simulated data sets in `.Rdata` format.

``` r
### Initialization of the parameters
sizeEffect <- log(3) #The size effect
N.effectif <- 200 #The sample size
N.stop <- 10000 #The numbert of simulated data sets
path <- "/home/Simulations/" # Work directory
Size.Effect <- "log3" #The folder names for the size effect
```

Note that to simulate the data sets of the simple scenario, only the `sim.base` function needs to be modified with the variables and coefficients described in Supplementary Table A2 of the article (see tables at the bottom of the page).


The second file `2_estimate_param_ATE.R` provides scripts to obtain the theoretical values of the marginal effects by the mean of the unadjusted estimations for 1,000,000 simulated data sets. As in the previous script, only the parameters at the beginning need to be specified to obtain the outcomes for different scenarios.


## Estimating the treatment effects according to the learners

We considered several models and algorithms (learners) to fit the outcome model. All analyses were performed using R version 4.3.0, using the `caret` package with a tuning grid of length equal to 20. Below an overview of the learners used:

- **Lasso logistic regression.** L1 regularization allows for the selection of the predictors. To establish a flexible model, we considered all possible interactions between the treatment arm A and the covariates X. Additionally, we used B-splines for the continuous covariates. The `glmnet` package was used. The penalization of the L1 norm was the only tuning parameter.

- **Elasticnet logistic regression.** This approach mirrors the logistic regression mentioned earlier but incorporates both the L1 and the L2 regularizations.

- **Neural network.** We chose one hidden layer, which represents one of the most common network architectures. Its size constitutes the tuning parameter. The `nnet` package was used.

- **Support vector machine.** To relax the linear assumption, we opted for the radial basis function kernel. The `svmRadial` function of the `kernlab` package was used. It requires two tuning parameters: the cost penalty of miss-classification and the flexibility of the classification.

- **Super learner.** We also tested a super learner with the ensemble of the previous ML techniques. Super learner consisted in a weighting average of the learner-specific predictions by using a weighted linear predictor. In alignment with our previous choices, we estimated the weights by maximizing the average AUC through a 20-fold cross-validation. We used the `SuperLearner` package.

The third script `3_complex_sim.R` applies all the previously described techniques to the simulated data sets. The script creates a new folder named `/results_complex_nXX_ateYY`. The output is a `.csv file` containing a list of estimates for each technique and for each simulated data set (as shown below).


To use this script on the simple scenario, the model specifications for both penalized methods need to be adjusted accordingly. This modification is needed in the function for the Superlearner, for the tuning parameters and inside the bootstrapping process to predict the estimates. See R code below with an example for obtaining the tuning parameters for the Elasticnet model in both the complex and the simple scenario.

``` r
### Elasticnet for the complex scenario
elasticnet.param <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) +
 bs(x5, df = 3) + 	x6 + x7 + bs(x8, df = 3) + x9 + bs(x10, df = 3) + bs(x11, df = 3) +
 x12 + bs(x14, df = 3) + bs(x15, df = 3) +	bs(x18, df = 3) + x19 + x20 + bs(x21, df = 3)),
 data = base.train, 	method = 'glmnet',  tuneLength = 20, metric = "ROC", trControl = control,
 family = "binomial", penalty.factor = c(0, rep(1, 78)))

### Elasticnet for the simple scenario
elasticnet.param <- train(outcome ~ ttt * (bs(x1, df = 3) + bs(x2, df = 3) + bs(x3, df = 3) +
 x4 + x5 + x6), data = base.train, method = 'glmnet',  tuneLength = 20, metric = "ROC",
 trControl = control, family = "binomial", penalty.factor = c(0,rep(1,24)))
``` 		   

This R script has a long processing time, so parallel processing is highly recommended. The code was written to allow for stopping and resuming execution from where it was left off.


## Estimating the performances of the learners


The fourth R script `4_read_sims.R` contains the functions necessary to compute the performance criteria for each of the learners. It creates a data frame with the learners in the rows and the performances criteria in the columns. Finally, it saves the table as a `.tex` file.



## Supplementary tables


| **Variable** | **Role in the study**  | **Distribution**  |
| :- | :--- | :-------- |
$X_1$   |  Continuous   covariate   |  $N (0,1)$     | 
$X_2$   |  Continuous    covariate  |       $N (0,1)$   | 
$X_3$    |  Continuous  covariate   |       $N (0,1)$    | 
$X_4$   |  Continuous   covariate  |      $1 \{ N (0,1)< -0.67  \}$  (i.e., prevalence $\sim$ 25\%)      | 
$X_5$   | Binary   covariate   |        $1 \{ N (0,1)< 0.67  \}$  (i.e., prevalence $\sim$ 50\%)  | 
$X_6$   |  Binary   covariate   |      ${1}  \{ N (0,1)< 0.67  \}$  (i.e., prevalence $\sim$ 75\%)   | 
$A$  | Binary treatment arm    |     ${1}  \{ N (0,1)> 0  \}$  (i.e., a 1:1 randomized clinical trial)  | 
$Y$  | Binary outcome  |  $B ( n, p =  logistic \; (  \beta_0 + \beta_1  \times X_1 + \beta_1  \times X_2 + \beta_1  \times X_4 + \beta_1  \times X_5 +   \beta_2 \times A  )   )$   |
: A2-Model used for the simulations in the simple scenario
