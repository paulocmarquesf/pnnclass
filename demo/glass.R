# Forensic glass dataset

library(pnnclass)
library(microbenchmark)
library(ggplot2)
library(extrafont)
library(latex2exp)

#       min       lq     mean   median      uq      max neval
#  8.805169 9.010958 9.825337 9.151921 9.68024 16.38385   100

# microbenchmark({
model <- pnnclass(type ~ ., data = Glass_trn) # , distance = "manhattan")

pred <- predict(model, newdata = Glass_tst)
#})

print(mean(pred$y_hat != Glass_tst$type), digits = 4)

table(pred$y_hat, Glass_tst$type, dnn = c("Predicted", "Observed"))

###

pca <- prcomp(scale(Glass_tst, center = TRUE, scale = TRUE))
cumsum(pca$sdev) / sum(pca$sdev)
pca$rotation[, 1:2]
db <- as.data.frame(pca$x[, 1:2])
db$y <- Glass_tst$type
db$y_hat <- pred$y_hat
db$hit <- ifelse(db$y == db$y_hat, "Yes", "No")

theme_set(theme_bw())

ggplot(db, aes(x = PC1, y = PC2, color = hit, shape = factor(y))) +
    geom_point(size = 1, show.legend = FALSE) +
    scale_colour_manual(values = c("red", "blue")) +
    labs(x = "Score on PC1",
         y = "Score on PC2",
         title = "Forensic Glass testing sample",
         color = "Correctly classified") +
    theme(text = element_text(family = "Times New Roman", size = 10))
