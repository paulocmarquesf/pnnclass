pnnclass <- function(formula, data, distance = "euclidean", beta_max = 10) {
    mf <- model.frame(formula, data)
    X_trn <- as.matrix(mf[-1])
    y_trn <- mf[[1]] # model.response(mf) is a named vector
  
    dt <- .data_topology(X_trn, y_trn, distance)
    
    n <- nrow(X_trn)
    beta_hat <- numeric(n - 1)
    log_p <- numeric(n - 1)
    for (r in 1:(n - 1)) {
        opt <- optim(par = 1,
                     fn = function(beta_r) - .log_p_r(beta_r, dt$A[r], dt$C[r, ], dt$L),
                     lower = 0,
                     upper = beta_max,
                     method = "L-BFGS-B")
        beta_hat[r] <- opt$par
        log_p[r] <- -opt$value
    }
    
    cv <- .loocv(dt, beta_hat)

    model <- list(predictors = attr(attr(mf, "terms"), "term.labels"),
                  dt = dt, beta_hat = beta_hat, log_p = log_p, Q = cv$Q,
                  loo_error = cv$error, kappa_hat = cv$kappa_hat)
    class(model) <- "pnnclass"
    model
}

predict.pnnclass <- function(model, newdata) {
    X_tst <- as.matrix(newdata[model$predictors])
    pred <- .predict_test(model, X_tst)
    pred$y_hat <- apply(pred$prob, 1, function(row) which.max(row) - 1)
    pred
}
