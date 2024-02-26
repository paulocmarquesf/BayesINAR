# BayesINAR

Bayesian predictive models for time series of counts.

> Paulo C. Marques F., Helton Graziadei and Hedibert F. Lopes.
"_Bayesian generalizations of the integer-valued autoregressive model_."
Journal of Applied Statistics, v. 1, p. 1-21, 2020. https://doi.org/10.1080/02664763.2020.1812544

```r
devtools::install_github("paulocmarquesf/BayesINAR")

library(BayesINAR)

y <- Pittsburgh$Area_51

model <- inar(y)
summary(model)
predict(model, h = 2)

model <- adinar(y)
summary(model)
predict(model, h = 2)

model <- dpinar(y)
summary(model)
predict(model, h = 2)

cross_validate(y, model = "inar")
```
