#' @export logisticEM
#' @import data.table
logisticEM <- function(y, id, x, beta, lambda, k = length(lambda), epsilon = 1e-06, maxit = 100) {
    if (ncol(beta) != k)
        stop("Column number of initial beta does not match class number!")
    if (length(lambda) != k)
        stop("Length of initial lambda does not match class number!")
    n <- length(y)
    p <- ncol(x)
    n_id <- length(unique(id))
    all.loglik <- -Inf
    iter <- 0

    while (T) {
        xbeta <- data.matrix(x) %*% beta
        temp <- y * xbeta - log1pexp(xbeta)
        temp <- as.data.table(temp)
        temp <- temp[, lapply(.SD, sum), id]
        logf <- as.matrix(temp)[, -1, drop = F]
        rid <- as.matrix(temp)[, 1]

        if (any(lambda <= 0))
            stop("lambda cannot be equal to or smaller than 0!")
        # if (any(lambda>=1)) stop('lambda cannot be equal to or larger than 1!') if (sum(lambda)!=0)
        # stop('lambda must sum up to 1!')
        if (requireNamespace("matrixStats", quietly = TRUE)) {
            RM <- matrixStats::rowMaxs(logf)
        } else {
            RM <- apply(logf, 1, max)
        }
        logf2 <- logf - RM
        temp <- t(lambda * t(exp(logf2)))
        loglik <- sum(log(rowSums(temp)) + RM)
        all.loglik <- c(all.loglik, loglik)
        improve <- diff(all.loglik)
        if (improve[length(improve)] < epsilon | iter >= maxit) break

        posteriorz <- temp/rowSums(temp)
        lambda.old <- lambda
        lambda <- colMeans(posteriorz)

        beta.old <- beta
        w <- match(id, rid)
        updatebeta <- lapply(1:k, function(j) glm(y ~ . - 1, data = data.frame(as.matrix(x)), weights = posteriorz[w,
            j], family = binomial()))
        beta <- sapply(updatebeta, coef)
        iter <- iter + 1
    }
    se.beta <- sapply(updatebeta, function(x) summary(x)$coef[, 2])
    z.beta <- beta/se.beta
    temp <- abs(z.beta)
    p.beta <- 2 - 2 * pnorm(temp)
    rownames(posteriorz) <- rid
    list(lambda = lambda, beta = beta, posteriorz = posteriorz, all.loglik = all.loglik[-1], se.beta = se.beta,
        p.beta = p.beta, y = y, id = id, x = x)
}
log1pexp <- function(x) {
    y <- x
    y[x <= -37] <- exp(x[x <= -37])
    y[x > -37 & x <= 18] <- log1p(exp(x[x > -37 & x <= 18]))
    y[x > 18 & x <= 33.3] <- x[x > 18 & x <= 33.3] + exp(-x[x > 18 & x <= 33.3])
    y
}
