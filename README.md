R codes used in the article titled: G-computation for increasing performances of randomized clinical trials with binary response
================

## Description

In the article (DOI. XX), we investigated by simulations the performances of G-computation according various machine learning algorithms in the context of randomized clinical trials with a 1:1 allocation ratio, a superiority design and binary outcomes. 

## Data generation

Two distinct scenarios were simulated. In the provided code we consider the complex scenario where the estimation of the outcome model was not straightforward. In the article, we considered different marginal treatment effects and sample sizes. For all the scenarios 10,000 data sets were generated.

[DAGcomplexe.pdf](https://github.com/user-attachments/files/17360824/DAGcomplexe.pdf)

| **Variable** | **Role in the Study** | **Distribution** |
|--------------|-----------------------|------------------|
| X₁           | Continuous covariate  | 𝓝(0, 1)          |
| X₂           | Continuous covariate  | 𝓝(β₀ + β₁ * X₁, 1)  |
| X₃           | Continuous covariate  | 𝓝(β₀ - β₁ * X₁ - β₂ * X₂, 1) |
| X₄           | Continuous covariate  | 𝓝(0, 1)          |
| X₅           | Binary covariate      | 𝟙{𝓝(0, 1) > 0.67} (i.e., prevalence ~ 25%)  |
| X₆           | Binary covariate      | 𝟙{𝓝(β₀ - β₁ * X₄, 1) > -0.40} (i.e., prevalence ~ 50%) |
| X₇           | Continuous covariate  | 𝓝(β₀ - β₁ * X₅, 1) |
| X₈           | Binary covariate      | 𝟙{𝓝(β₀ + β₁ * X₆, 1) > -0.80} (i.e., prevalence ~ 75%) |
| X₉           | Continuous covariate  | 𝓝(β₀ + β₁ * X₇, 1) |
| X₁₀          | Continuous covariate  | 𝓝(0, 1)          |
| X₁₁          | Binary covariate      | 𝟙{𝓝(β₀ + β₁ * X₈, 1) > 0.84} (i.e., prevalence ~ 25%) |
| X₁₂          | Continuous covariate  | 𝓝(β₀ - β₁ * X₁₁ - β₂ * X₁₀, 1) |
| X₁₃          | Continuous covariate  | 𝓝(β₀ - β₁ * X₁₁, 1) |
| X₁₄          | Continuous covariate  | 𝓝(0, 1)          |
| X₁₅          | Binary covariate      | 𝟙{𝓝(0, 1) > 0.67} (i.e., prevalence ~ 25%) |
| X₁₆          | Binary covariate      | 𝟙{𝓝(0, 1) > 0.67} (i.e., prevalence ~ 25%) |
| X₁₇          | Continuous covariate  | 𝓝(0, 1)          |
| A            | Binary treatment arm   | 𝟙{𝓝(0, 1) > 0} (i.e., a 1:1 randomized clinical trial) |
| Y            | Binary outcome        | 𝓑(n, p = logistic(β₂ + β₃ * 𝟙{X₂ > -0.44} - β₃ * X₃ + (β₃ / 2) * X₃² + β₃ * X₅ + β₃ * X₆ + β₃ * X₉ + (β₃ / 2) * X₁₀² - β₃ * X₁₂ - β₃ * (X₁₃ > -0.55) + β₃ * X₁₄ + β₃ * X₁₅ + (β₃ / 2) * A * X₁₄ + β₄ * A) |
* 𝟙{a} = 1 if the condition a is true and 0 otherwise; 𝓝(μ, σ) represents a Gaussian distribution with mean at μ and standard deviation at σ; 𝓑(n, p) represents a Binomial distribution with a size n and probability of success p. The regression coefficients were: β₀ = -0.4, β₁ = log(2), β₂ = -2, β₃ = log(2), and β₄ = log(3); log(1.5); log(0.9729) to obtain mOR values of 1.9; 1.3; 1.00 respectively.



## Machine learning techniques

We considered several models and algorithms (learners) to fit the outcome model. All analyses were performed using R version 4.3.0, using the caret package with a tuning grid of length equal to 20. Below an overview of the learners used:

• Lasso logistic regression. L1 regularization allows for the selection of the predictors. To establish a flexible model, we considered all possible interactions between the treatment arm A and the covariates X. Additionally, we used B-splines for the continuous covariates. The glmnet package was used. The penalization of the L1 norm was the only tuning parameter.
• Elasticnet logistic regression. This approach mirrors the logistic regression mentioned earlier but incorporates both the L1 and the L2 regularizations.
• Neural network. We chose one hidden layer, which represents one of the most common network architectures. Its size constitutes the tuning parameter. The nnet package was used.
• Support vector machine. To relax the linear assumption, we opted for the radial basis function kernel. The svmRadial function of the kernlab package was used. It requires two tuning parameters: the cost penalty of miss-classification and the flexibility of the classification.
• Super learner. We also tested a super learner with the ensemble of the previous ML techniques. Super learner consisted in a weighting average of the learner-specific predictions by using a weighted linear predictor. In alignment with our previous choices,17 we estimated the weights by maximizing the average AUC through a 20-fold cross-validation. We used the SuperLearner package.

``` r
install.packages("test")
# Output test
#> Likelihood ratio test=70.5213 on 3 df, p=0
#> n=1000, number of events=369
```

This code is a specific case of the complex scenario for n=200 and a mOR=1.9