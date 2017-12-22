#' DCmods.
#'
#' This package is a collection of various statistical models that can prove to
#' be useful for marketing and social research on choice behavior data. The
#' models included are the changepoint model, consideration set model, and
#' latent stage model.
#'
#' For details on each model, please refer to the package's vignette, "Better
#' Understanding DCmods," as well as the help pages of each function.
#'
#' @name DCmods
#' @docType package
#' @author Elizabeth Bruch, Fred Feinberg, Cameron Hollingshead, Yue Wang,
#'   Yujing Zhou
#'
#' @import chngpt
#' @importFrom changepoint cpt.mean cpt.var cpt.meanvar penalty_decision logLik
#' @import data.table
#' @importFrom stats pchisq binomial coef glm model.matrix optim pnorm rnorm
#' @importFrom zoo coredata
NULL
