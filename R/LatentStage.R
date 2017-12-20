#' fit discrete heterogeneity model with multiple stages
#'
#' Discrete heterogeneity model applies to scenario each subject is assigned a latent
#' class, and fixed parameters are used across subjects in the same class. Multiple stages
#' simply mean that more than one responses are allowed in fitting this model. Typically
#' one response is nested with another one, but it doesn't have to be this case to fit the
#' model.
#'
#' Arbitrary number of covariates, repsonses and subject ids can be givin in the argument
#' list, however, their number must be the same, which equal to the number of stages in
#' this model. And of course, their dimension must match, which means row(X1) = length(y1)
#' = length(id1).
#'
#' @param nclass numbers of class for all subjects, this determines how many groups of
#' different parameters will be obtained
#' @param X1 covariate matrix of the first stage, similar for X2, X3, ...
#' @param y1 responsed of the first stage, similar for y2, y3, ...
#' @param id1 subject id of the first stage, similar for id2, id3, ...
#' @return The output is a list object containing both arguments and results. \code{lambda}
#' is the estimate of class proportions, which sum up to 1. \code{beta} contains estimates, SEs
#' and p-values for all linear parameters. \code{posteriorz} lists the probability of each subject
#' belonging to a specific group, the estimated class id is determined by which class maximizes
#' the probability. \code{all.loglik} lists log-likelihood in each iteration. \code{y} is a list
#' object of responses of all stages. \code{id} is a list object of subject ids of all stages.
#' \code{x} is a list object of covariates of all stages. \code{AIC} is the AIC for current
#' model. \code{BIC} is the BIC for current model. These two can be used for class number
#' selection. \code{runtime} is the whole elapsed time to fit the model.
#' @examples
#' data(threestage)
#' attach(threestage)
#' mod <- LatentStage(5, X1=stage1[, 4:7],
#'        y1=stage1$Y1, id1=stage1$Person,
#'        X2=stage2[, 4:7],
#'        y2=stage2$Y2, id2=stage2$Person,
#'        X3=stage3[, 4:7],
#'        y3=stage3$Y3, id3=stage3$Person)
#'
#' data(dating)
#' attach(dating)
#' nonmiss <- !is.na(wrote)
#' mod <- LatentStage(3, y1 = browsed, y2 = wrote[nonmiss],
#'                    id1 = respid, id2 = respid[nonmiss],
#'                    X1 = agedif, X2 = agedif[nonmiss])
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
    # printit <- c('beta', 'lambda', 'BIC', 'runtime')
    # print(out[printit])
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
