#' @title Max test for the joint association between multiple traits
#' and environmental variables  (incl. nominal variables) using weighted averaging methods (CWM/SNC regression)
#'
#' @description
#' For the test on multi-trait multi-environment association, two test statistics are available, one which accounts
#' for correlations among the traits and among the environmental variables (score_test == TRUE) and one
#' which does not. The former is related to double constrained correspondence analysis (ter Braak et al. 2018) and the latter to RLQ (Dray et al. 2014).
#' created using \code{\link{make_obj_for_traitenv}}. If E is a TE_obj, the arguments L and T are ignored.
#' @param type 'joint' testing of the association between a set of traits and a set environmental variables (default) or
#' 'onebyone' 1-1 tests for each pair of trait and environmental variable
#' @param score_test logical. TRUE for the score test statistic, which accounts for correlations (default) and
#'  FALSE without accounting for correlations among traits and among environmental variables
#' @param fast logical, whether the fast version is used or the slow version (default TRUE).
#' For multi-trait-environment data, fast is set to TRUE.
#' @inheritParams CWMSNC_regressions
#' @details
#' Nomimal traits and environmental variables (factors) are expanded into dummy variables
#' if E and/or T are data frames. Set \code{E} or \code{T} to a data frame with a single factor
#' to obtain test of association with the factor.
#' For multi-trait, multi-environment data and score_test==TRUE, the trait and environment data
#' are internally transformed to normalized singular vectors. If score_test== FALSE, the variables are standardized
#' internally to weighted mean 0 and variance 1.
#' The permutation tests is carried
#' out by the max test. This test takes the maximum of the p-values of the CWM-E regression and
#' the SNC-T regression. Both p-values are obtained by permutation.
#' In the former, the rows of E are permuted, in the latter the rows of T are permuted.
#' The slow version recomputes the CWM and SNC and the associated regression from scratch for each permutation
#' and is available only with single numberic trait and environmental variable.
#' The fast version precomputes the CWM and SNC and computes the test statistic by a matrix product in which
#' the permuted traits or environmental variables are rescaled to weigthed mean 0 and variance 1 using the
#' precomputed weights (determined by the weighing option).
#' @return  A named list, among which,
#' \item{p_values}{five p-values: permutational, site-based and species-based and its maximum, and the parametric ones}
#' \item{obs}{values of the test statistic for sites (first row) and species (second row)}
#' \item{sim.row}{values of the test statistic for the nrepet data in which the rows of E are permuted}
#' \item{sim.col}{values of the test statistic for the nrepet data in which the rows of T are permuted}
#' \item{nrepet}{number of permutations}
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#'
#' ter Braak, C.J.F., Šmilauer, P. & Dray, S. (2018) Algorithms and biplots for double constrained correspondence analysis. Environmental and Ecological Statistics, 25, 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x
#'
#' Dray, S., Choler, P., Dolédec, S., Peres-Neto, P.R., Thuiller, W., Pavoine, S. & ter Braak, C.J.F. (2014) Combining the fourth-corner and the RLQ methods for assessing trait responses to
#' environmental variation. Ecology, 95, 14-21. http://dx.doi.org/10.1890/13-0196.1
#' @example demo/TutorialWA_Aravo_Multi.r
#' @export WA_p_max
WA_p_max <- function(E,L=0,T=0, nrepet = 19,  weighing = "N2", type= 'joint', score_test = TRUE, cutoff = 0, fast = TRUE) {
  if (!type %in% c("joint","onebyone")) warning("type in WA-p_max is joint or onebyone")
  # local function: make a  design matrix (with for factors all levels included) for name x in data dat
  dummy_coding <- function(x, dat){model.matrix(as.formula(paste("~ 0", paste(x,  collapse= "+"), sep= "+")), data = dat)}
  #end local function
  if(class(E)[1]!="TE_obj"){
    # E and T might be factors or data frames with factors
    if (is.factor(E) ||is.factor(T)) {
      if (is.factor(E)) E <- as.data.frame(E)
      if (is.factor(T)) T <- as.data.frame(T)
    }
    if (is.data.frame(E)) {
      if (type == "joint"){
        E <-  model.matrix(as.formula(paste("~ 0", paste(names(E),  collapse= "+"), sep= "+")), data = E)
      } else { # onebyone (make sure each factor has all levels included in the design matrix)
        # combine all design matrices
        E<-Reduce(cbind,lapply(names(E), function(x){dummy_coding(x,E)}))
      }
    }
    if (is.data.frame(T)) {
      if (type == "joint"){
        T <-  model.matrix(as.formula(paste("~ 0", paste(names(T),  collapse= "+"), sep= "+")), data = T)
      } else { # onebyone
        T<-Reduce(cbind,lapply(names(T), function(x){dummy_coding(x,T)}))
      }
    }
    obj <- make_obj_for_traitenv(E, L, T, cutoff)
  } else obj <- E
  if (!fast && ncol(obj$E) + ncol(obj$T) != 2 ) {
    warning("option fast set to TRUE, because there are more than one trait or environmetnal variable")
    fast <- TRUE
  }
  if (!fast){
    # following the description with test stat = 1/p_value of regressions performed in full
    result <- WA_p_max0(obj, nrepet = nrepet,  weighing = weighing)
  } else {
    # short cut, which is even more general (allows multi-trait multi-envi)
    if(ncol(obj$E)==1 && ncol(obj$T == 1)) score_test <- FALSE
    result <- WA_p_max1(obj, nrepet = nrepet,  weighing = weighing, type = type, wsvd = score_test)
  }
  return(result)
}
