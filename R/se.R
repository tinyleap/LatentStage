expit <- function(x) {
    1/(1 + exp(-x))
}
#' calculate outer-version standard errors of linear parameters and class proportions
#'
#' Calculate outer-version standard errors for estimate of linear parameters
#' and class proportions
#'
#' @param out a list object obtained directly from \code{\link{LatentStage}}
#' @return a list object with two attributes. \code{se.beta} is the outer-version SE of
#' linear parameters. \code{se.lambda} is the outer-version SE of class proportions. Note
#' in the \code{\link{LatentStage}} function, the latter is not provided.
#' @export
se_outer <- function(out) {
    xbeta <- data.matrix(out$x) %*% out$beta
    temp <- out$y * xbeta - log1pexp(xbeta)
    temp <- as.data.table(temp)
    temp <- temp[, lapply(.SD, sum), out$id]
    logf <- as.matrix(temp)[, -1, drop = F]
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        RM <- matrixStats::rowMaxs(logf)
    } else {
        RM <- apply(logf, 1, max)
    }
    logf2 <- logf - RM
    f <- exp(logf2)
    L <- c(f %*% matrix(out$lambda))
    f <- f/L

    nclass <- length(out$lambda)
    dlogL <- NULL
    for (k in 1:nclass) {
        temp <- (out$y - expit(xbeta[, k])) * out$x
        temp <- as.data.table(as.matrix(temp))
        temp <- temp[, lapply(.SD, sum), out$id]
        dlogf_dbeta <- as.matrix(temp)[, -1]
        dlogL <- cbind(dlogL, out$lambda[k] * f[, k] * dlogf_dbeta)
    }
    jacobi <- Matrix::bdiag(diag(dim(dlogL)[2]), diag(out$lambda) - tcrossprod(out$lambda))  # pi to lambda
    jacobi2 <- Matrix::bdiag(diag(dim(dlogL)[2]), rbind(diag(nclass - 1), -1))  # lambda to lin.ind
    jacobi <- as.matrix(jacobi)
    jacobi2 <- as.matrix(jacobi2)
    nb <- 1:dim(dlogL)[2]
    np <- dim(dlogL)[2] + 1:dim(f)[2]
    dlogL <- cbind(dlogL, f)
    dlogL <- dlogL %*% jacobi
    dlogL <- dlogL %*% jacobi2

    N <- length(unique(out$id))  # what on earth is N?
    B <- N/(N - 1) * t(dlogL) %*% dlogL
    Sigma <- jacobi2 %*% solve(B) %*% t(jacobi2)  # back to lin.dep
    se <- sqrt(diag(Sigma))
    list(se.beta = matrix(se[nb], ncol = nclass, dimnames = dimnames(out$beta)), se.lambda = se[np])
}
se_standard <- function(out){
    xbeta <- out$x %*% out$beta
    expitxbeta <-expit(xbeta)
    nperson <- max(out$id)
    nclass <- length(out$lambda)
    nvar <- dim(out$beta)[1]
    dlogf_dbeta <- array(0, c(nperson, nclass, nvar))
    d2logf_dbeta2 <- array(0, c(nperson, nclass, nvar^2))
    xtx <- apply(out$x, 1, tcrossprod)
    xtx <- t(xtx)
    for (k in 1:nclass){
        temp <- expitxbeta[, k] * (1 - expitxbeta[, k]) * xtx
        temp <- as.data.table(temp)
        temp <- temp[, lapply(.SD, sum), out$id][, -1]
        temp <- as.matrix(temp)
        d2logf_dbeta2[, k,] <- (-temp)
        temp <- (out$y - expitxbeta[, k]) * out$x
        temp <- as.data.table(as.matrix(temp))
        temp <- temp[, lapply(.SD, sum), out$id][, -1]
        temp <- as.matrix(temp)
        dlogf_dbeta[, k,] <- temp
    }
}
