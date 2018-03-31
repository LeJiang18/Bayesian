# Bayesian Logistic Regression

**Project Goal**

Implement Bayesian Logistic Regression and use it to analyze the Diabetes dataset.

**Idea**

Given the prior distribution of parameter, and Polya-Gamma distribution, developed the marginal posterior distribution of the parameter of interest using Gibbs Sampler. We chose zero and MLE as the initial value for the hyper-parameter, and compared the ESS, then determined the sample size and burn in size. We did sensitivity analysis to see the difference as a result of different starting point, and the outcome an informative and a flat distribution can lead to. Lastly, in the variable selection section, we dropped most insignificant variables and variables that are highly correlated.
