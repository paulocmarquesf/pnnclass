# pnnclass

Probabilistic nearest neighbors classification

> Bruno Fava, Paulo C. Marques F., and Hedibert F. Lopes. 2024. "Probabilistic Nearest Neighbors Classification" Entropy 26, no. 1: 39. https://doi.org/10.3390/e26010039

```r
devtools::install_github("paulocmarquesf/pnnclass")

library(pnnclass)
```

```r
model <- pnnclass(y ~ x1 + x2, data = Ripley_trn)
pred <- predict(model, newdata = Ripley_tst)
round(mean(pred$y_hat != Ripley_tst$y), 4)
```

```
0.0840
```

```r
model <- pnnclass(type ~ ., data = MASS::Pima.tr)
pred <- predict(model, newdata = MASS::Pima.te)
round(mean(pred$y_hat != MASS::Pima.te$type), 4)
```

```
0.2199
```

```r
model <- pnnclass(type ~ ., data = Glass_trn)
pred <- predict(model, newdata = Glass_tst)
round(mean(pred$y_hat != Glass_tst$type), 4)
```

```
0.2804
```
