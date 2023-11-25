pnnclass <- function(formula, data, distance = "euclidean", beta_max = 10) {
    mfrm <- model.frame(formula, data)
    
    if (!all(apply(mfrm[-1], 2, class) == "numeric")) stop("All predictors must be numeric.")
    
    X_trn <- as.matrix(mfrm[-1])
    y_trn <- as.integer(as.factor(mfrm[[1]])) - 1
    y_levels <- levels(as.factor(mfrm[[1]]))

    dt <- .data_topology(X_trn, y_trn, distance)
    
    n <- nrow(X_trn)
    beta_hat <- numeric()
    for (r in 1:(n - 1)) {
        opt <- optim(par = 1, fn = function(beta_r) - .log_p_r(beta_r, dt$A[r], dt$C[r, ], dt$L),
                     lower = 0, upper = beta_max, method = "L-BFGS-B")
        if (opt$par == 0) break
        beta_hat <- c(beta_hat, opt$par)
    }
    
    cv <- .loocv(dt, beta_hat)
    
    model <- list(predictors = attr(attr(mfrm, "terms"), "term.labels"),
                  y_levels = y_levels,
                  beta_max = beta_max, dt = dt, beta_hat = beta_hat,
                  Q = cv$Q, loo_error = cv$error, k_hat = cv$k_hat)
    
    class(model) <- "pnnclass"
    
    model
}

predict.pnnclass <- function(model, newdata) {
    X_tst <- as.matrix(newdata[model$predictors])
    pred <- .predict_test(model, X_tst)
    pred$y_hat <- model$y_levels[apply(pred$prob, 1, function(row) which.max(row))]
    pred
}
