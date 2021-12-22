# Forensic glass dataset

library(pnnclass)

model <- pnnclass(type ~ ., data = Glass_trn)

model$kappa_hat
model$loo_error

pred <- predict(model, newdata = Glass_tst)

mean(pred$y_hat != Glass_tst$type)
