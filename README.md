```bibtex
@article{marquesf2022,
  author = {Paulo C. {Marques F.} and Helton Graziadei and Hedibert F. Lopes},
  title = {Bayesian generalizations of the integer-valued autoregressive model},
  journal = {Journal of Applied Statistics},
  volume = {49},
  number = {2},
  pages = {336--356},
  year = {2022},
  publisher = {Taylor \& Francis},
  doi = {https://doi.org/10.1080/02664763.2020.1812544},
}
```

## Bayesian generalizations of the integer-valued autoregressive model

> Paulo C. Marques F., Helton Graziadei and Hedibert F. Lopes

> https://doi.org/10.1080/02664763.2020.1812544

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