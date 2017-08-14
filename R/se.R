expit <- function(x) {
    1/(1 + exp(-x))
}
se_outer <- function(out) {
    xbeta <- out$xb %*% out$betab
    temp <- out$yb * xbeta - log1pexp(xbeta)
    temp <- as.data.table(as.matrix(temp))
    temp <- temp[, lapply(.SD, sum), out$idb]
    logf <- as.matrix(temp)[, -1]
    logf2 <- logf - rowMaxs(logf)  # It is such a genius of me!
    f <- exp(logf2)
    L <- c(f %*% matrix(out$lambda))
    f <- f/L

    nclass <- length(out$lambda)
    dlogL <- NULL
    for (k in 1:nclass) {
        temp <- (out$yb - expit(xbeta[, k])) * out$xb
        temp <- as.data.table(as.matrix(temp))
        temp <- temp[, lapply(.SD, sum), out$idb]
        dlogf_dbeta <- as.matrix(temp)[, -1]
        dlogL <- cbind(dlogL, out$lambda[k] * f[, k] * dlogf_dbeta)
    }
    jacobi <- bdiag(diag(dim(dlogL)[2]), diag(out$lambda) - tcrossprod(out$lambda))  # pi to lambda
    jacobi2 <- bdiag(diag(dim(dlogL)[2]), rbind(diag(nclass - 1), -1))  # lambda to lin.ind
    nb <- 1:dim(dlogL)[2]
    np <- dim(dlogL)[2] + 1:dim(f)[2]
    dlogL <- cbind(dlogL, f)
    dlogL <- dlogL %*% jacobi
    dlogL <- dlogL %*% jacobi2

    N <- max(out$idb)  # what on earth is N?
    B <- N/(N - 1) * t(dlogL) %*% dlogL
    Sigma <- jacobi2 %*% solve(B) %*% t(jacobi2)  # back to lin.dep
    se <- sqrt(diag(Sigma))
    list(se.beta = matrix(se[nb], ncol = nclass, dimnames = dimnames(out$betab)), se.lambda = se[np])
}
# se_standard <- function(out){
#     xbeta <- out$xb %*% out$betab
#     expitxbeta <-expit(xbeta)
#     nperson <- max(out$idb)
#     nclass <- length(out$lambda)
#     nvar <- dim(out$betab)[1]
#     dlogf_dbeta <- array(0, c(nperson, nclass, nvar))
#     d2logf_dbeta2 <- array(0, c(nperson, nclass, nvar^2))
#     xtx <- apply(out$xb, 1, tcrossprod)
#     xtx <- t(xtx)
#     for (k in 1:nclass){
#         temp <- expitxbeta[, k] * (1 - expitxbeta[, k]) * xtx
#         temp <- as.data.table(temp)
#         temp <- temp[, lapply(.SD, sum), out$idb][, -1]
#         temp <- as.matrix(temp)
#         d2logf_dbeta2[, k,] <- (-temp)
#         temp <- (out$yb - expitxbeta[, k]) * out$xb
#         temp <- as.data.table(as.matrix(temp))
#         temp <- temp[, lapply(.SD, sum), out$idb][, -1]
#         temp <- as.matrix(temp) dlogf_dbeta[, k,] <- temp
#     }
# }
