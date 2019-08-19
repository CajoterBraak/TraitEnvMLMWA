#' @title Permutational max test of trait-environment assocation for MLM models
#'
#' @description
#' \code{MLM3_p_max} performs the permutational max test of trait-environment interaction for MLM models,
#' starting from a fitted MLM3 object.
#' @param  test_stat choice of test statistic; 'Wald' (default) or 'LRT'.
#' @param  how_to_permute a list for two \code{\link[permute]{how}} calls, the first for site permutation, the second for species permutation.
#' @param print integer; print progress every \code{print} iterations
#' @inheritParams Bootstrap_test_prmtrc
#' @details The code assumes that the interaction parameter is the last fixed parameter in summary(MLM3).
#' The data used in the max test is extracted using \code{\link{dat4MLM2TE_obj}}. This generates an object of class \code{TE_obj} (see
#' \code{\link{make_obj_for_traitenv}}). \code{\link{dat4MLM2TE_obj}} is limitted for use with a single trait and single environmental variable,
#' and so is therefore \code{MLM3_p_max}.
#' In the model-based permutation tests, either the trait values or the environmental values in the interaction term T*E of the model are permuted to yield
#' a species- and site-level test, respectively.
#' Main effects for the permuted trait and environmental variable are added to ensure
#' that the interaction after permutation has a corresponding main effect.
#' For further details, see Appendices A4 and A5 of ter Braak (2019).
#' @return  A named list, among which,
#' \item{p_values}{four p-values: one parametric p-value (Wald test) and three permutational p-values: site-based and species-based and the maximum of these two values}
#' \item{obs}{values of the test statistic for sites (first row) and species (second row)}
#' \item{sim.row}{values of the test statistic for the nrepet data in which the rows of E are permuted}
#' \item{sim.col}{values of the test statistic for the nrepet data in which the rows of T are permuted}
#' \item{nrepet}{number of permutations}
#' @seealso \code{\link{expand4glmm}}.
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @examples
#' \dontrun{
#' #use a precomputed MLM3 model for the Revisit data
#' data("MLM3")
#' ## or compute an MLM3 model from the data
#' # data("Revisit")
#' # formula.MLM3 <- y ~ poly(env,2) + poly(trait,2) +
#'   env : trait  + (1 + env|species) + (1 + trait| site)
#' # MLM3 <- glmmTMB(formula.MLM3, family = betabinomial,  data=Revisit)
#' summary(MLM3)

#' res_perm <- MLM3_p_max(MLM3, test_stat = "Wald", nrepet = nrepet, Binomial_total = 100)
#' names(res_perm)
#' round(res_perm$p_values,3)
#' }
#' @export

MLM3_p_max <- function(MLM3, nrepet = 19, Binomial_total = 0, test_stat = "Wald", how_to_permute=list(sites=how(),species=how()), print =1, nAGQ=0) {
  #n_sites <- nrow(obj$L); n_species <- ncol(obj$L)

  options1 <- setoptions4pmax_test(verbose = TRUE,
            FUN_test_statistics =  list(init= init_MLM_pmax, test= MLM_Wald_bte),
            with.MLM= TRUE, nrepet = nrepet,
            perm.mat.spe = 0 , perm.mat.sit = 0, print = print,
            nmax = 6, how_to_permute=how_to_permute, filek = 1)

  if (test_stat == "Wald") {with.LRT<-FALSE; LRT <- FALSE} else {with.LRT<-TRUE; LRT <- TRUE; }

  if (class(MLM3)[1] =="glmmTMB") {library = "glmmTMB";B <- summary(MLM3)$coefficients$'cond'  } else {library = "lme4";B <- summary(MLM3)$coefficients}
  options2 <- setoptions4MLM( formula = formula(MLM3),
                              family = family(MLM3), K = Binomial_total,
                              model.names= c("MLM3"), estimation = FALSE,
                              library = library, output = "short", nAGQ =nAGQ,
                              with.LRT= with.LRT, LRT.only = LRT, verbose = FALSE)
  options3 <- c(options1,options2[-length(options2)])
  obj <- dat4MLM2TE_obj(model.frame(MLM3)) # TE_obj
  if (test_stat == "Wald")  obj$p_parametric <- B[nrow(B),4] else obj$p_parametric <- NA
  result <- PermutationTest_r_c_max(obj, options = options3)
  result$p_values <- c(p_parametric = result$p_parametric[1], result$p_values)
  return(result)
}

