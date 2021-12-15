# devtools::install_github("paulocmarquesf/BayesINAR")

library(BayesINAR)

data(Pittsburgh)
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

# cross_validate(y, model = "inar")
