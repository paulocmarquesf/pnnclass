```bibtex
@article{fava2024,
  title = {Probabilistic Nearest Neighbors Classification},
  author = {Bruno Fava and Paulo C. {Marques F.} and Hedibert F. Lopes},
  journal = {Entropy},
  volume = {26},
  year = {2024},
  number = {1},
  issn = {1099-4300},
  doi = {https://doi.org/10.3390/e26010039},
  url = {https://www.mdpi.com/1099-4300/26/1/39}
}
```

# Probabilistic nearest neighbors classification

> Bruno Fava, Paulo C. Marques F., and Hedibert F. Lopes

> https://doi.org/10.3390/e26010039

**Abstract.** Analysis of the currently established Bayesian nearest neighbors classification model points to a connection between the computation of its normalizing constant and issues of NPcompleteness. An alternative predictive model constructed by aggregating the predictive distributions of simpler nonlocal models is proposed, and analytic expressions for the normalizing constants of these nonlocal models are derived, ensuring polynomial time computation without approximations. Experiments with synthetic and real datasets showcase the predictive performance of the proposed predictive model.

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