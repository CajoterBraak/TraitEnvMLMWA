# Appendix to ter Braak 2019
#New robust weighted averaging- and model-based methods for assessing trait-environment relationships.
MLM_Wald_bte <- function(obj, options,  ...){
  #  function for use in PermutationTest_r_c_max
  # @param obj object of class TE_obj (made by function make_obj_for_traitenv(E, L, T,cutoff))
  # @formula model formula using names created by expand4glmm(obj, K = K)[use verbose = TRUE to see these names]
  # @param K is binomial total; if K[1]==0 data are not binomial or presence / absence
  #                      K can be a scalar or vector of n*m (number of sites * number of species)
  # @param family glm family; use character form
  # @param library use "glmmTMB" or else (e.g. "lme4")
  # @param verbose If TRUE, prints the names for use in formula.
  # @param estimation if TRUE, the interaction coefficient is given, else testing is carried out via the Wald test
  # @value test statistic: 1/pWald_b_te
  result <- with(options,
                 {
                   if (estimation){
                     col_b_te = 1 # column position of the interaction coefficient to select it
                     invert_p_value = FALSE
                   } else col_b_te = 4   # column position of the p-value to select it
                   dat <- expand4glmm(obj, K = K)
                   #names(dat)[c(5:8)] <- c("trait", "env","trait0", "env0")
                   #if (verbose) print(str(dat))
                   if (library == "glmmTMB"){
                     suppressWarnings( MLM <- try(glmmTMB(formula, data=dat, family= family)))
                     if (!class(MLM)[1]=="try-error"){
                       B <- summary(MLM)$coefficients$'cond'}
                     else {B <- matrix(NA, nrow =2, ncol = 4); MLM = 0}
                   } else if  (library == "lme4") {
                     suppressWarnings(MLM <- try(glmer(formula, data=dat, family= family,  nAGQ=nAGQ, control = glmerControl(calc.derivs=F))))
                     if (!class(MLM)[1]=="try-error"){
                       B <- summary(MLM)$coefficients
                     } else {B <- matrix(NA, nrow =2, ncol = 4); MLM = 0}
                   } else {print(paste("library", library ,"not implemented"))}
                   # BEWARE the interaction coefficient is likely the last one.... adapt if not
                   p_val_Wald <- B[nrow(B),col_b_te]
                   if (verbose){
                     print(summary(MLM))
                     if (estimation) {names(p_val_Wald) <- "b_te"} else names(p_val_Wald) <- "p_val_Wald"
                     print(p_val_Wald)
                   }
                   if (is.character(family) ) testnam <- paste(library, family, "MLMWald", sep = ".") else testnam <- paste(library, "MLMWald", sep = ".")
                   if (invert_p_value & !estimation) test_stat <- 1/p_val_Wald else test_stat <- p_val_Wald
                   #test_stat <- matrix(test_stat, nrow = 2, ncol = 1, dimnames =list(c("rows","cols"), testnam) )
                   if (with.LRT){
                     formulaNULL <- update(formula(MLM),  ~ . - trait:env )
                     if (class(MLM)[1] =="glmmTMB") {
                       suppressWarnings(MLM0 <- try(update(MLM, formula = formulaNULL)))
                     } else {
                       suppressWarnings(MLM0 <- try(update(MLM, formulaNULL,  nAGQ=nAGQ, control = glmerControl(calc.derivs=F))))
                     }
                     if (class(MLM0)[1]=="try-error") {
                       if (class(MLM)[1] =="glmmTMB") suppressWarnings(MLM0 <- try(glmmTMB(formulaNULL, data=dat, family= family(MLM)))) else
                         suppressWarnings(MLM0 <- try(glmer(formula, data=dat, family= family,  nAGQ=nAGQ, control = glmerControl(calc.derivs=F))))
                       if (class(MLM0)[1]=="try-error") {
                         p_val_LRT <- NA
                       } else {
                         p_val_LRT <- pchisq(2*( logLik(MLM) - logLik(MLM0)), df =1, lower.tail = FALSE)
                       }
                     } else p_val_LRT <- pchisq(2*( logLik(MLM) - logLik(MLM0)), df =1, lower.tail = FALSE)
                     if (invert_p_value) test_stat2 <- 1/p_val_LRT else test_stat2 <- p_val_LRT
                     names(test_stat2) <- "pval LRT"
                     if (is.character(family) ) testnam2 <- paste(library, family, "MLM.LRT", sep = ".") else testnam2 <- paste(library, "MLM.LRT", sep = ".")

                     #test_stat2 <- matrix(test_stat2, nrow = 2, ncol = 1, dimnames =list(c("rows","cols"), testnam2) )
                     # print(test_stat2)
                     if (LRT.only) {
                       test_stat <- test_stat2
                     } else {
                       test_stat <- c(test_stat,test_stat2)
                       names(test_stat) <- c("p_val_Wald", "p_val_LRT")
                     }
                   }
                   #print(test_stat)
                   if (estimation) {
                     test_stat <- c(test_stat, B[nrow(B),2]) # standard error
                     names(test_stat) <- c("b_te","se_b_te")
                   } #else rownames(test_stat)<- c("teststat.site","teststat.species")
                   if (verbose) print(test_stat)
                   if (output =="short"){
                     result <- test_stat
                   } else {
                     result <- list(test_stat = test_stat, model = MLM )
                   }
                   return(result)
                 })
  return(result)
}


init_MLM_pmax <- function(obj, options){list(obj=obj, options=options)}
