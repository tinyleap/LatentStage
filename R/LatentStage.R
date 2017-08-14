LatentStage <- function(nclass, ...) {
    tt <- proc.time()
    dots <- list(...)
    dots <- dots[order(names(dots))]
    y <- dots[grep("y\\d", names(dots))]
    obsid <- dots[grep("id\\d", names(dots))]
    X <- dots[grep("X\\d", names(dots))]

    yall <- do.call(c, y)
    obsidall <- do.call(c, obsid)
    X <- lapply(X, function(x) model.matrix(~., x))
    names(X) <- NULL
    Xall <- do.call(bdiag, X)

    Xall <- Xall[order(obsidall), ]
    yall <- yall[order(obsidall)]
    obsidall <- obsidall[order(obsidall)]
    out <- logisticEM(y = yall, id = obsidall, x = Xall, beta = matrix(rnorm(dim(Xall)[2] * nclass), ncol = nclass),
                      lambda = rep(1, nclass)/nclass)
    out$runtime <- proc.time() - tt
    out
}
getmaxperm <- function(est_cls, cls) {
    if(max(est_cls) != max(cls)) stop('estimated class number is wrong!')
    pm <- permutations(max(cls), max(cls))
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
