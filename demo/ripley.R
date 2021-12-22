# Ripley synthetic dataset

library(pnnclass)

model <- pnnclass(y ~ x1 + x2, data = Ripley_trn)

model$kappa_hat
model$loo_error

pred <- predict(model, newdata = Ripley_tst)

mean(pred$y_hat != Ripley_tst$y)
