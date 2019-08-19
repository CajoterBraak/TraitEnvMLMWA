#' @title Parametric bootstrap test of trait-environment assocation using the MLM models
#'
#' @description
#' \code{Bootstrap_test_prmtrc} performs a parametric bootstrap test on trait-environment interaction,
#' starting from a fitted MLM3 object (or any other MLM object with the trait:env term as last fixed parameter).
#' @param  MLM3 the fitted MLM3 object, created by glmer (lme4) or glmmTMB.
#' @param  test_stat choice of test statistic; 'Wald' (default), 'LRT' or 'both'. The default is quicker.
#' @param  nrepet number of bootstraps
#' @param  Binomial_total scalar, 0 for count-like data and
#'  the binomial total for logit models (1 for presence-absence).
#' @param nAGQ  integer scalar (default 0), used only for an object created by glmer
#' @details The code assumes that the parameter for trait-environment interaction (trait:env) is the last fixed parameter in summary(MLM3).
#' First, the formula of the null model is created by deleting the trait:env term from the formula of object MLM3 (the non-null model).
#' Second, the null model (MLM0) is fitted. Three, sampling is from this null model.
#' For each simulated data set, the model is refitted using the formula of the MLM3 object.
#' The code works therefore also for MLM1 and MLM2, although their use is not recommended.
#' See also Box A2 in Appendix A4 and Appendix A1.
#' @return  A named list,
#' \item{p_values}{the parametric and bootstrap p-values}
#' \item{MLM0}{the fitted null model}
#' \item{obs}{value of the test statistic(s)}
#' \item{nrepet}{number of bootstraps}
#' \item{sim.boot}{values of the test statistic(s) for the nrepet bootstrapped data}
#' \item{test_stat}{the chosen test statistic(s)}
#' @seealso \code{\link{expand4glmm}}.
#' @import lme4 glmmTMB stats
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @examples
#' \dontrun{
#' #use a precomputed MLM3 model, e.g. from the Revisit data
#' data("MLM3")
#' ## or compute the MLM3 model from the data
#' # data("Revisit")
#' # formula.MLM3 <- y ~ poly(env,2) + poly(trait,2) +
#'  env : trait  + (1 + env|species) + (1 + trait| site)
#' # MLM3 <- glmmTMB(formula.MLM3, family = betabinomial,  data=Revisit)
#' summary(MLM3)
#' res_boot <- Bootstrap_test_prmtrc(MLM3, test_stat = "Wald", nrepet = nrepet, Binomial_total = 100)
#' names(res_boot)
#' round(res_boot$p_values,3)
#' }
#' @export
Bootstrap_test_prmtrc <- function(MLM3, test_stat = "Wald", nrepet = 19, Binomial_total = 0, nAGQ = 0) {
  # Bootstrap test based on test_stat (Wald or LRT or both) according to Box A2 in Appendix A3.
  # MLM3 is a fitted glmmTMB or lme4 objects
  # Binomial_total  = 0 for Poisson/Negative binomial counts and Binomial total for (beta-)binomial
  # BEWARE the interaction coefficient is likely the last one, in the code nrow(B), .... adapt if not

  # Step 1: arguments test_stat (Wald or LRT) and nrepet
  # Step 2: fit the null model
  formulaNULL <- update(formula(MLM3),  ~ . - trait:env )
  dat <- model.frame(MLM3)
  if (class(MLM3)[1] =="glmmTMB") {
  #  suppressWarnings(MLM0 <- glmmTMB(formulaNULL, data=dat, family= family(MLM3)))
    suppressWarnings(MLM0 <- update(MLM3, formulaNULL))
  } else {
    suppressWarnings(MLM0 <- update(MLM3, formulaNULL,  nAGQ=nAGQ, control = glmerControl(calc.derivs=F)))
  }
  # Step 3: get F0: the test statistic from the data
  if (test_stat %in% c("Wald", "both")) {
    if (class(MLM3)[1] == "glmmTMB"){
      B <- summary(MLM3)$coefficients$'cond'
    } else { # lme4
      B <- summary(MLM3)$coefficients
    }
    # BEWARE the interaction coefficient is likely the last one.... adapt if not
    obs1 <- 1/B[nrow(B),4] # 1/p-value
    names(obs1) <- "Wald"
  }
  if (test_stat %in% c("LRT", "both")) {
    obs2 <-  logLik(MLM3) - logLik(MLM0)
    names(obs2) <- "LRT"
  }
  if (test_stat=="Wald") obs <- obs1 else if (test_stat=="LRT") obs <- obs2 else obs<- c(obs1,obs2)
  if (test_stat == "both") id_LRT <-2 else id_LRT <- 1
  obj <- dat4MLM2TE_obj(dat) # TE_obj
  n_sites <- nrow(obj$L); n_species <- ncol(obj$L)

  sim.boot <- matrix(0, nrow = nrepet, ncol = length(obs))
  colnames(sim.boot) <- paste("boot",names(obs), sep ="")
  for (i in 1:nrepet){
    # step 4 and 5: simulated from the null model and calculate the test statistic
    Ysim <- simulate(MLM0)
    L <- matrix(Ysim[[1]], nrow = n_sites,ncol = n_species, dimnames = list(sites= rownames(obj$L),species=colnames(obj$L)))
    obj.sim <- make_obj_for_traitenv(obj$E,L,obj$T)
    dat.sim <- expand4glmm(obj.sim, K = Binomial_total)
    if (class(MLM3)[1] == "glmmTMB"){
      suppressWarnings(MLM3.sim <- try(glmmTMB(formula(MLM3), data=dat.sim, family= family(MLM3))))
      if (!class(MLM3.sim)[1]=="try-error"){
        B <- summary(MLM3.sim)$coefficients$'cond'}
      else {B <- matrix(NA, nrow =2, ncol = 4); MLM3.sim = 0}
    } else if (class(MLM3)[1] == "glmerMod") {
      suppressWarnings(MLM3.sim <- try(glmer(formula(MLM3), data=dat.sim, family= family(MLM3),  nAGQ=nAGQ, control = glmerControl(calc.derivs=F))))
      if (!class(MLM3.sim)[1]=="try-error"){
        B <- summary(MLM3.sim)$coefficients
      } else {B <- matrix(NA, nrow =2, ncol = 4); MLM3.sim = 0}
    }
    if (test_stat %in% c("Wald","both")){
      sim.boot[i,1] <- 1/B[nrow(B),4]
    }
    if (test_stat %in% c("LRT","both")) { # LRT: refit null model
      if (class(MLM3.sim)[1] =="glmmTMB") {
        suppressWarnings(MLM0.sim <- try(glmmTMB(formulaNULL, data=dat.sim, family= family(MLM3))))
      } else {
        suppressWarnings(MLM0.sim <- try(update(MLM3.sim, formulaNULL,  nAGQ=nAGQ, control = glmerControl(calc.derivs=F))))
      }
      if (class(MLM0)[1]=="try-error") {
        sim.boot[i,id_LRT] <- NA
      } else sim.boot[i,id_LRT] <- logLik(MLM3.sim) - logLik(MLM0.sim)
    }
  }
  # step 6: 	Compute the Monte Carlo significance level,
  #           i.e. compute the number of values F 0, F 1, F 2, ... , FK
  #           that is greater than or equal to F0 (this number is thus at least 1),
  #           and divide by K + 1 (with K = nrepet).
  if (length(obs)==1){
    isna.r <-  sum(is.na(sim.boot))
    p_boot <- (sum(abs(sim.boot) >= abs(obs), na.rm = TRUE) + 1)  / (nrepet- isna.r  + 1)
  } else {
    obs.mat <- matrix(rep(abs(obs),each=nrepet), nrow= nrepet, ncol=length(obs))
    isna.r <-  colSums(is.na(sim.boot))
    p_boot <- (colSums(abs(sim.boot) >= obs.mat, na.rm=TRUE) + 1)/ (nrepet - isna.r + 1)
  }
  if (test_stat %in% c("Wald","both")) p_val1 <- 1/obs[1]
  if (test_stat %in% c("LRT","both"))  p_val2 <- pchisq(2*obs[id_LRT], df =1, lower.tail = FALSE)
  if (test_stat=="Wald") p_val <- p_val1 else if (test_stat=="LRT") p_val <- p_val2 else p_val<- c(p_val1,p_val2)
  if (test_stat == "both") {
    result <- cbind(p_prmtrc= p_val, p_boot = p_boot)

  } else result <- c(p_prmtrc= p_val, p_boot = p_boot)
  attr(result, "nrepet")<- nrepet
  return(list(p_values=result, nrepet = nrepet, sim.boot=sim.boot, MLM0 = MLM0, obs = obs, test_stat = test_stat))
}
