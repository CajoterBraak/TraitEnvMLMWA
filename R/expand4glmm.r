#' @title Utility function expand4glmm
#'
#' @description
#' \code{expand4glmm} inflates a \code{TE_obj} object (see \code{\link{make_obj_for_traitenv}}), which consists of
#' trait, environment and abundance data, into a format for glm, glmer and glmmTMB.
#' @param  obj an object of class TE_obj, usually the output of function \code{\link{make_obj_for_traitenv}}.
#' @param  K scalar, 0 for count-like data and
#'  the binomial total for data that you like to analyse using logit models (1 for presence-absence).
#' @details With single trait and single environment variable data (ncol(obj$T)== ncol(obj$E)==1), the names
#' of the trait and environmental variables are changed to trait and env (to allow identical MLM formulas for different data and variables ).
#' In the multi-variable case,
#' the original names are kept and all interactions are added.
#' \code{expand4glmm} is used repeatedly in the model-based permutational max test
#'  using function \code{\link{MLM3_p_max}}. In this function, \code{obj} is extended
#'  with list elements \code{trait0} and \code{env0}
#'  to allow for the detail needed in the model-based permutational max test. If used in this way, the
#'  internal logical \code{inPermut_r_c} is set to \code{TRUE}.
#' @return  A data frame. For single trat and enviornmental variable data:
#' \item{y}{the response; the abuncance values in obj$L }
#' \item{site}{a factor indicating the site corresponding to each value of y}
#' \item{species}{a factor indicating the species corresponding to each value of y }
#' \item{obs}{integer 1 to N, the length of y}
#' \item{trait}{trait}
#' \item{env}{environmental variable}
#' \item{trait.env}{the product}
#' In the multi-trait environment case:
#' all traits and environmental variables and all pairwise interactions
#' @examples
#' \dontrun{
#' data("aravo",  package = "ade4")
#' Y <-aravo$spe
#' trait <- scale(aravo$traits$SLA)
#' env <- scale(aravo$env$Snow)
#' obj <- make_obj_for_traitenv(env,Y, trait, cutoff=0)
#' K  = 0 # 0 for non-binomial models
#' dat <- expand4glmm(obj, K = K)
#' names(dat)
#' library(glmmTMB)
#' formula.MLM3 <- y ~ poly(env,2) + poly(trait,2) +
#'  env : trait  + (1 + env|species) + (1 + trait| site)
#' MLM3 <- glmmTMB(formula.MLM3, family = betabinomial,  data=dat)
#' summary(MLM3)
#' plot_MLM3(MLM3)
#' }
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278 )
#'
#' ter Braak, C.J.F., Peres-Neto, P. & Dray, S. (2017) A critical issue in
#' model-based inference for studying trait-based community assembly and a solution.
#' PeerJ, 5, e2885. https://doi.org/10.7717/peerj.2885
#'
#' @export expand4glmm
expand4glmm <- function(obj, K= 0){
  # K  = 0 for Poisson/Negative binomial counts and Binomial total for (beta-)binomial
  # adapted from Jamil et al 2013
  if ('trait0' %in% names(obj)) inPermut_r_c = TRUE else inPermut_r_c = FALSE
  with(obj, {
    sitespec <- expand.grid(rownames(L),colnames(L))
    site <-sitespec[,1]; species<-sitespec[,2]
    y <- as.vector(as.matrix(L))
    Evec <-  E[site,, drop = FALSE];  rownames(Evec)= NULL
    Tvec <-  T[species,, drop = FALSE]; rownames(Tvec)= 1:nrow(Tvec)
    if (inPermut_r_c){
      Evec0 <-  env0[site,, drop = FALSE];  rownames(Evec0)= NULL
      Tvec0 <-  trait0[species,, drop = FALSE]; rownames(Tvec0)= 1:nrow(Tvec0)
    } else {Evec0 =NULL; Tvec0 = NULL}
    if (is.null(colnames(E)))  {
      colnames(Evec) = paste("env", 1:ncol(Evec), sep = "")
    }
    if (is.null(colnames(T))) {
      colnames(Tvec) = paste("trait", 1:ncol(Tvec), sep = "")
    }
    ET <- Rten2_with_names(Evec,Tvec)
    if (inPermut_r_c){ # ncol(ET)==1: single trait/env
      XYZ <- data.frame(y,site,species, obs = 1:length(y),Tvec, Evec, Tvec0,  Evec0, ET)
      if (ncol(ET)==1) names(XYZ)[c(5:8)] <- c("trait", "env","trait0", "env0") else {
        # not one trait ,env
        print(names(XYZ))
        stop("more than one trait or env detected in expand4glmm which is currently not supported in MLM3_p_max")
        # adapt formulas for more than one trait-env
      }
    } else {
      XYZ <- data.frame(y,site,species, obs = 1:length(y),Tvec, Evec, ET)
      if (length(names(XYZ)==7))names(XYZ)[5:6] <- c("trait", "env") # length(names(XYZ)==6): single trait/env
    }
    if (!is.factor(site)) XYZ$site <- factor(XYZ$site)
    if (!is.factor(species)) XYZ$species <- factor(XYZ$species)
    if (K[1]>0) XYZ$y = cbind(y, K-y)
    return(XYZ) # XYZ notation of Jamil et al 2013
  })
}
Rten2_with_names <- function(X1,X2) {
  one.1 <- matrix(1,1,ncol(X1))
  one.2 <- matrix(1,1,ncol(X2))
  res  <- kronecker(X1,one.2)*kronecker(one.1,X2) #  column of X2 runs fastest
  colnames(res) <- kronecker(colnames(X1),colnames(X2),function(a,b){ paste(b, a, sep = ":")} )
  rownames(res) <- rownames(X1)
  return(res)
}
