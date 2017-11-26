#' @export
LatentStage <- function(nclass, ...) {
    tt <- proc.time()
    dots <- list(...)
    dots <- dots[order(names(dots))]
    y <- dots[grep("y\\d", names(dots))]
    obsid <- dots[grep("id\\d", names(dots))]
    X <- dots[grep("X\\d", names(dots))]
    if (length(y)!=length(obsid)) stop("Stage of y and id do not match!")
    if (length(y)!=length(X)) stop("Stage of y and X do not match!")
    X <- lapply(X, function(x) model.matrix(~., data.frame(x)))
    for (i in 1:length(y)) {
        if (!all(y[[i]] %in% c(0, 1))) stop('y should only take 0 or 1!')
        if (length(y[[i]])!=length(obsid[[i]])) stop('Length of y and id do not match!')
        if (length(y[[i]])!=dim(X[[i]])[1]) stop('Dimension of y and X do not match!')
    }

    yall <- unlist(y)
    obsidall <- unlist(obsid)
    Xall <- Matrix::bdiag(X)
    covname <- sapply(X, colnames)
    covname <- paste('stage', col(covname), covname)


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
        beta <- beta[, od, drop=F]
        posteriorz <- posteriorz[, od, drop=F]
        se.beta <- se.beta[, od, drop=F]
        p.beta <- p.beta[, od, drop=F]
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
    temp <- sapply(out[c('beta', 'se.beta', 'p.beta')], function(x)c(t(x)))
    rownames(temp) <- apply(expand.grid(colnames(out$beta), rownames(out$beta)), 1, function(x)Reduce(paste, x))
    out$beta <- temp
    out$se.beta <- NULL
    out$p.beta <- NULL
    printit <- c('beta', 'lambda', 'BIC', 'runtime')
    print(out[printit])
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
