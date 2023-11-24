# Pima indian diabetes dataset

library(pnnclass)
library(pROC)
library(ggplot2)
library(extrafont)
library(latex2exp)

#       min       lq     mean   median       uq      max neval
#  52.67951 57.28314 59.64538 59.04403 60.04796 114.9671   100

# microbenchmark::microbenchmark({
model <- pnnclass(type ~ ., data = MASS::Pima.tr)

pred <- predict(model, newdata = MASS::Pima.te)
# })

print(mean(pred$y_hat != MASS::Pima.te$type), digits = 4)

###

roc_pnn <- roc(MASS::Pima.te$type, pred$prob[, 2])

par(family = "Times New Roman", cex = 0.85)

plot(roc_pnn, col = "blue", grid = TRUE,
     xlab = "FPR (1 - Specificity)",
     ylab = "TPR (Sensitivity)",
     main = "ROC Pima",
     print.thres = "best",
     print.thres.pattern = "Threshold = %.3f\nSpecificity = %.3f\nSensitivity = %.3f",
     print.auc = TRUE, print.auc.adj = c(0.5, 5), print.auc.col = "black",
     legacy.axes = TRUE, asp = FALSE, las = 1)

auc(roc_pnn)

idx_best <- with(roc_pnn, which.max(sensitivities + specificities))
roc_pnn$thresholds[idx_best]
roc_pnn$specificities[idx_best]
roc_pnn$sensitivities[idx_best]
