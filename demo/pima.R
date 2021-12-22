# Pima indian diabetes dataset

library(pnnclass)

model <- pnnclass(type ~ ., data = Pima_trn)

model$kappa_hat
model$loo_error

pred <- predict(model, newdata = Pima_tst)

mean(pred$y_hat != Pima_tst$type)
