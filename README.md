R scripts for data simulation in de Keizer et al. (submitted paper, 2024)
================


## Description

We aimed to investigate by simulations the performances of G-computation according to various machine learning algorithms in the context of randomized clinical trials with a 1:1 allocation ratio, a superiority design and binary outcomes. 


## Simulating the 10,000 datasets

This repository contains all the R scripts for the complex scenario with a sample size of 200 and a marginal odds ratio (mOR) of 1.9. The code is optimized for the Linux environment but is also executable on Windows.  At the beginning of each R script, a list of parameters is available for simulating the other complex scenarios shown in the article. The first file `1_sim_data.R` contains the scripts to simulate data sets with the following complex scenario. 

<p align="center">
<img src="https://github.com/user-attachments/assets/699a9861-48c8-4bd3-8d9b-88f01d25a9fe" alt="image" width="700" height="auto">
</p>

Only the first lines need to be modified to run the code with different effect and sample sizes. When executed, this script will create a folder named `/complex_nXX_ateYY`, where `XX` is the sample size (`N.effectif`) and `YY` the effect size (`SizeEffect`), containing all the simulated data sets in `.Rdata` format.

``` r
### Initialization of the parameters
SizeEffect <- "log(3)" #The size effect
N.effectif <- 200 #The sample size
N.stop <- 10000 #The number of simulated data sets
path <- "/home/Simulations/" # Work directory
```

Note that to simulate the data sets of the simple scenario, only the `sim.base` function needs to be modified with the variables and coefficients described in Supplementary Table A2 of the article (see tables at the bottom of the page).


The second file `2_estimate_param_ATE.R` provides scripts to obtain the theoretical values of the marginal effects by the mean of the unadjusted estimations for 1,000,000 simulated data sets. As in the previous script, only the parameters at the beginning need to be specified to obtain the outcomes for different scenarios.


## Estimating the treatment effects according to the learners

We considered several models and algorithms (learners) to fit the outcome model. All analyses were performed using R version 4.3.0, using the `caret` package with a tuning grid of length equal to 20. Below an overview of the learners used:

- **Lasso logistic regression.** L1 regularization allows for the selection of the predictors. To establish a flexible model, we considered all possible interactions between the treatment arm A and the covariates X. Additionally, we used B-splines for the continuous covariates. The `glmnet` package was used. The penalization of the L1 norm was the only tuning parameter.

- **Elasticnet logistic regression.** This approach mirrors the logistic regression mentioned earlier but incorporates both the L1 and the L2 regularizations.

- **Neural network.** We chose one hidden layer, which represents one of the most common network architectures. Its size constitutes the tuning parameter. The `nnet` package was used.

- **Support vector machine.** To relax the linear assumption, we opted for the radial basis function kernel. The `svmRadial` function of the `kernlab` package was used. It requires two tuning parameters: the cost penalty of miss-classification and the flexibility of the classification.

- **Super learner.** We also tested a super learner with the ensemble of the previous machine learning techniques. Super learner consisted in a weighting average of the learner-specific predictions by using a weighted linear predictor. In alignment with our previous choices, we estimated the weights by maximizing the average area under the receiver operating characteristic curve through a 20-fold cross-validation. We used the `SuperLearner` package.

The third script `3_complex_sim.R` applies all the previously described techniques to the simulated data sets. The script creates a new folder named `/results_complex_nXX_ateYY`. The output is a `.csv file` containing a list of estimates for each technique and for each simulated data set (as shown below).


To use this script on the simple scenario, the model specifications for both penalized methods need to be adjusted accordingly. This modification is needed in the function for the Super learner, for the tuning parameters and inside the bootstrapping process to predict the estimates. See R code below with an example for obtaining the tuning parameters for the Elasticnet model in both the complex and the simple scenario.

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


### Table A1 : Models used for simulations in the complex scenario

| **Variable** | **Role in the study** | **Distribution** |
|--------------|-----------------------|------------------|
| $X_1$        | Continuous covariate  | $\mathcal{N}(0, 1)$          |
| $X_2$        | Continuous covariate  | $\mathcal{N}(\beta_0 + \beta_1 \times X_1, 1)$               |
| $X_3$        | Continuous covariate  | $\mathcal{N}(\beta_0 - \beta_1 \times X_1 - \beta_2 \times X_2, 1)$               |
| $X_4$        | Continuous covariate  | $\mathcal{N}(0, 1)$               |
| $X_5$        | Binary covariate      | $\mathbf{1}\{\mathcal{N}(0, 1) > 0.67\}$ (i.e., prevalence $\sim$ 25\%)              |
| $X_6$        | Binary covariate      | $\mathbf{1}\{\mathcal{N}(\beta_0 - \beta_1 \times X_4, 1) > -0.40\}$ (i.e., prevalence $\sim$ 50\%)         |
| $X_7$        | Continuous covariate  | $\mathcal{N}(\beta_0 - \beta_1 \times X_5, 1)$         |
| $X_8$        | Binary covariate      | $\mathbf{1}\{\mathcal{N}(\beta_0 + \beta_1 \times X_6, 1) > -0.80\}$ (i.e., prevalence $\sim$ 75\%)    |
| $X_9$        | Continuous covariate  | $\mathcal{N}(\beta_0 + \beta_1 \times X_7, 1)$  |
| $X_{10}$     | Continuous covariate  | $\mathcal{N}(0, 1)$ |
| $X_{11}$     | Binary covariate      | $\mathbf{1}\{\mathcal{N}(\beta_0 + \beta_1 \times X_8, 1) > 0.84\}$ (i.e., prevalence $\sim$ 25\%)  |
| $X_{12}$     | Continuous covariate  | $\mathcal{N}(\beta_0 - \beta_1 \times X_{11} - \beta_2 \times X_{10}, 1)$ |
| $X_{13}$     | Continuous covariate  | $\mathcal{N}(\beta_0 - \beta_1 \times X_{11}, 1)$ |
| $X_{14}$     | Continuous covariate  | $\mathcal{N}(0, 1)$ |
| $X_{15}$     | Binary covariate      | $\mathbf{1}\{\mathcal{N}(0, 1) > 0.67\}$ (i.e., prevalence $\sim$ 25\%)  |
| $X_{16}$     | Binary covariate      | $\mathbf{1}\{\mathcal{N}(0, 1) > 0.67\}$ (i.e., prevalence $\sim$ 25\%)  |
| $X_{17}$     | Continuous covariate  | $\mathcal{N}(0, 1)$    | 
| $A$          | Binary treatment arm  | $\mathbf{1}\{\mathcal{N}(0, 1) > 0\}$ (i.e., a 1:1 randomized clinical trial)   |
| $Y$          | Binary outcome        | $\mathcal{B}\big(n, p = \text{logistic}\big(\beta_2 + \beta_3 \times \mathbf{1}\{X_{2} > -0.44\} - \beta_3 \times X_3 + (\beta_3 / 2) \times X_3^2$ |
|              |                       | $+ \beta_3 \times X_5 + \beta_3 \times X_6 + \beta_3 \times X_9 + (\beta_3 / 2) \times X_{10}^2 - \beta_3 \times X_{12}$ |
|              |                       | $- \beta_3 \times (X_{13} > -0.55) + \beta_3 \times X_{14} + \beta_3 \times X_{15} + (\beta_3 / 2) \times A \times X_{14} + \beta_4 \times A \big)\big)$ |


**Notes**:
$\mathbf{1}\{ a \} = 1$ if the condition $a$ is true and 0 otherwise; $\mathcal{N}(\mu, \sigma)$ represents a Gaussian distribution with mean at $\mu$ and standard deviation at $\sigma$; $\mathcal{B}(n, p)$ represents a Binomial distribution with size $n$ and probability of success $p$.
The regression coefficients were: $\beta_0 = -0.4$, $\beta_1 = \log(2)$, $\beta_2 = -2$, $\beta_3 = \log(2)$, and $\beta_4 = \log(3), \log(1.5), \log(0.9729)$ to obtain mOR values of 1.9, 1.3, and 1.0 respectively.




### Table A2 : Models used for simulations in the simple scenario

| **Variable** | **Role in the study** | **Distribution** |
|--------------|-----------------------|------------------|
| $X_1$        | Continuous covariate  | $\mathcal{N}(0, 1)$          |
| $X_2$        | Continuous covariate  | $\mathcal{N}(0, 1)$          |
| $X_3$        | Continuous covariate  | $\mathcal{N}(0, 1)$          |
| $X_4$        | Continuous covariate  | $\mathbf{1}\{\mathcal{N}(0, 1) < -0.67\}$ (i.e., prevalence $\sim$ 25\%) |
| $X_5$        | Binary covariate      | $\mathbf{1}\{\mathcal{N}(0, 1) < 0.67\}$ (i.e., prevalence $\sim$ 50\%) |
| $X_6$        | Binary covariate      | $\mathbf{1}\{\mathcal{N}(0, 1) < 0.67\}$ (i.e., prevalence $\sim$ 75\%) |
| $A$          | Binary treatment arm  | $\mathbf{1}\{\mathcal{N}(0, 1) > 0\}$ (i.e., a 1:1 randomized clinical trial) |
| $Y$          | Binary outcome        | $\mathcal{B}\big( n, p = \text{logistic} \big(  \beta_0 + \beta_1  \times X_1 + \beta_1  \times X_2 + \beta_1  \times X_4 + \beta_1  \times X_5 +   \beta_2 \times A    \big )  \big )$ |

**Notes**:
$\mathbf{1}\{  a \} = 1$ if the condition $a$ is true and 0 otherwise;  $\mathcal{N} (\mu,\sigma)$ represents a Gaussian distribution with mean $\mu$ and standard deviation $\sigma$; $\mathcal{B} (n,p)$ represents a Binomial distribution with a size $n$ and probability of success  $p$.  
The regression coefficients were: $\beta_0 = -3$,   $\beta_1 = \log(4)$  and $\beta_2 = \log(3), \log(1.5), \log(1)$ to obtain mOR values of 1.9, 1.3, 1.0 respectively.
For the simple scenario with reduced predictive performance the regression coefficients were: $\beta_0 = -0.6$,   $\beta_1 = 0.2$  and $\beta_2 = \log(3), \log(1.5), \log(1)$ to obtain mOR values of 3.0, 1.5, 1.0 respectively.









