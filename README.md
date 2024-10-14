R codes used in the article titled: G-computation for increasing performances of randomized clinical trials with binary response
================

## Description

In the article (DOI. XX), we investigated by simulations the performances of G-computation according various machine learning algorithms in the context of randomized clinical trials with a 1:1 allocation ratio, a superiority design and binary outcomes. 

## Data generation

Two distinct scenarios were simulated. In the provided code we consider the complex scenario where the estimation of the outcome model was not straightforward. In the article, we considered different marginal treatment effects and sample sizes. For all the scenarios 10,000 data sets were generated.

![dagcomplexe](https://github.com/user-attachments/assets/82e7b683-0814-45fd-8ca2-1b675d03dfc0)


![tablea1](https://github.com/user-attachments/assets/98ef0978-80c2-4038-aa8d-7d4bfc9c08f3)



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