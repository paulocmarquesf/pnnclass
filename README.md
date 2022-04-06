# pnnclass

Probabilistic nearest neighbors classification through aggregation of non-local models.

> Paulo C. Marques F. and Bruno Fava

```r
devtools::install_github("paulocmarquesf/pnnclass")

library(pnnclass)
```

```r
model <- pnnclass(y ~ x1 + x2, data = Ripley_trn)
pred <- predict(model, newdata = Ripley_tst)
mean(pred$y_hat != Ripley_tst$y)
```

```
0.084
```

```r
model <- pnnclass(type ~ ., data = Pima_trn)
pred <- predict(model, newdata = Pima_tst)
mean(pred$y_hat != Pima_tst$type)
```

```
0.2168675
```

```r
model <- pnnclass(type ~ ., data = Glass_trn, distance = "manhattan")
pred <- predict(model, newdata = Glass_tst)
mean(pred$y_hat != Glass_tst$type)
```

```
0.2429907
```
