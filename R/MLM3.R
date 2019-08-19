#' @title MLM3-fit to Revisit data using trait C:N ratio and environmental variable TMG
#'
#'
#' @description \code{MLM3} is a glmmTMB object consisting of the fit of the MLM3 model to the Revisit data using the
#' code
#'
#' \code{formula.MLM3 <- y ~ trait*env +(1+trait|site)+(1+env|species)}
#'
#' \code{MLM3 <- glmmTMB(formula.MLM3, family = betabinomial, data= Revisit)}
#' @seealso \code{\link{Revisit}}.
#' @references
#' ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)

#' @name MLM3
#' @docType data
NULL
