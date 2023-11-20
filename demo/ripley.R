# Ripley synthetic dataset

library(pnnclass)
library(pROC)
library(extrafont)
library(ggplot2)
library(latex2exp)

#       min       lq     mean   median       uq      max neval
#  118.8271 122.2817 126.8604 124.4563 125.9305 182.3407   100

# microbenchmark::microbenchmark({
model <- pnnclass(y ~ x1 + x2, data = Ripley_trn)

pred <- predict(model, newdata = Ripley_tst)
# })

print(mean(pred$y_hat != Ripley_tst$y), digits = 4)

###

roc_pnn <- roc(Ripley_tst$y, pred$prob[, 2])

# par(mfrow = c(1, 2))

par(family = "Times New Roman", cex = 0.85)

plot(roc_pnn, col = "blue", grid = TRUE,
     xlab = "FPR (1 - Specificity)",
     ylab = "TPR (Sensitivity)",
     main = "ROC Ripley",
     print.thres = "best",
     print.thres.pattern = "Threshold = %.3f\nSpecificity = %.3f\nSensitivity = %.3f",
     print.auc = TRUE, print.auc.adj = c(0.5, 5), print.auc.col = "black",
     legacy.axes = TRUE, asp = FALSE, las = 1)

auc(roc_pnn)

idx_best <- which.max(roc_pnn$sensitivities + roc_pnn$specificities)
roc_pnn$thresholds[idx_best]
roc_pnn$specificities[idx_best]
roc_pnn$sensitivities[idx_best]

###

grid <- expand.grid(x1 = seq(-1.5, 1.25, length.out = 100),
                    x2 = seq(-0.5, 1.5, length.out = 100))

pred <- predict(model, newdata = grid)

grid$pr = pred$prob[, 2]

theme_set(theme_bw())

ggplot() +
    scale_color_gradient(low = "cyan", high = "orange") +
    geom_point(data = grid, aes(x = x1, y = x2, color = pr),  show.legend = TRUE) +
    geom_point(data = Ripley_trn, aes(x = x1, y = x2, shape = factor(y)), size = 0.75, show.legend = FALSE) +
    labs(x  = TeX("$x_1$"), y = TeX("$x_2$"), title = "Ripley training sample and classifier probabilities") +
    # guides(shape = "none") +
    theme(text = element_text(family = "Times New Roman", size = 10),
          legend.position = "none")
 