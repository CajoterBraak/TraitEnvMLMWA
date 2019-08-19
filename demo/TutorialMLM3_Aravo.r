
# read and process data ---------------------------------------------------
data("aravo",  package = "ade4")
Y <-aravo$spe
trait <- scale(aravo$traits$SLA)
env <- scale(aravo$env$Snow)

obj <- make_obj_for_traitenv(env,Y, trait, cutoff=0)
str(obj)
K  = 0 # 0 = not binomial
dat <- expand4glmm(obj, K = K)
# note that the formula s below should only contain names available in dat, thus:
names(dat)

library(glmmTMB)

formula.MLM3 <- y ~ poly(env,2) + poly(trait,2) +  env : trait  + (1 + env|species) + (1 + trait| site)
MLM3 <- glmmTMB(formula.MLM3, family = "nbinom2",  data=dat)
summary(MLM3)

#Overdispersion parameter for nbinom2 family (): 1.4e+14

# switching to lme4 -------------------------------------------------------
library(lme4)
MLM3 <- glmer(formula.MLM3, family = poisson,  data=dat, nAGQ=0, control = glmerControl(calc.derivs=F))
summary(MLM3)

# could main effects be linear?
formula.MLM3.linear <- y ~ env * trait  + (1 + env|species) + (1 + trait| site)
MLM3.linear <- update(MLM3,formula.MLM3.linear)
anova(MLM3, MLM3.linear)


# Bootstrap confidence interval -------------------------------------------

getInteractionParameter <- function(MLM3){
  if (class(MLM3)[1] == "glmmTMB"){
    B <- summary(MLM3)$coefficients$'cond'
  } else if (class(MLM3)[1] == "glmerMod") {
    B <- summary(MLM3)$coefficients
  }
  # estimate with variance for boot.ci
  coef_te <- c(B[nrow(B),1],B[nrow(B),2]^2)
  return(coef_te)
}

aa <- bootMer(MLM3, FUN = getInteractionParameter, nsim = 1000)

library(boot)

b_te_boot <- aa$t[,1]
mean(b_te_boot, na.rm = TRUE)
sd(b_te_boot, na.rm = TRUE)
boot::boot.ci(aa, index = c(1,2), type = "stud" )



