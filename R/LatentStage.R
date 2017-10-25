#' @export LatentStage
LatentStage <- function(nclass, ...) {
    tt <- proc.time()
    dots <- list(...)
    dots <- dots[order(names(dots))]
    y <- dots[grep("y\\d", names(dots))]
    obsid <- dots[grep("id\\d", names(dots))]
    X <- dots[grep("X\\d", names(dots))]

    yall <- unlist(y)
    obsidall <- unlist(obsid)
    X <- lapply(X, function(x) model.matrix(~., data.frame(x)))
    covname <- c(sapply(X, colnames))
    Xall <- Matrix::bdiag(X)

    Xall <- Xall[order(obsidall), ]
    yall <- yall[order(obsidall)]
    obsidall <- obsidall[order(obsidall)]
    out <- logisticEM(y = yall, id = obsidall, x = Xall, beta = matrix(rnorm(dim(Xall)[2] * nclass), ncol = nclass),
                      lambda = rep(1, nclass)/nclass)
    od <- order(out$lambda)
    od <- rev(od)
    out <- within(out, {
        runtime <- proc.time() - tt
        lambda <- lambda[od]
        beta <- beta[, od]
        posteriorz <- posteriorz[, od]
        se.beta <- se.beta[, od]
        p.beta <- p.beta[, od]
        rownames(beta) <- covname
        rownames(se.beta) <- covname
        rownames(p.beta) <- covname
        names(lambda) <- paste('class', 1:nclass)
        colnames(beta) <- paste('class', 1:nclass)
        colnames(se.beta) <- paste('class', 1:nclass)
        colnames(p.beta) <- paste('class', 1:nclass)
        BIC <- -2*all.loglik[length(all.loglik)] + (length(beta) + nclass - 1)*log(length(unique(obsidall)))
        AIC <- -2*all.loglik[length(all.loglik)] + (length(beta) + nclass - 1)*2
    })
    out
}
getmaxperm <- function(est_cls, cls) {
    if(max(est_cls) != max(cls)) stop('estimated class number is wrong!')
    pm <- gtools::permutations(max(cls), max(cls))
    maxtr <- 0
    mt <- table(est_cls, cls)
    for (i in 1:dim(pm)[1]) {
        temp <- mt[pm[i, ], ]
        tr <- sum(diag(temp))
        if (tr > maxtr) {
            maxtr <- tr
            maxpm <- pm[i, ]
        }
    }
    maxpm
}
