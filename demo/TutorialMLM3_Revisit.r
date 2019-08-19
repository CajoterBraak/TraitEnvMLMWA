# get Revisit data ---------------------------------------------------
data("Revisit")

# parametric MLM analyses -------------------

library(glmmTMB)

## adapt dataframe Revisit for binomial analyses using  glm, glmmTMB or lme4
Revisit$y <- with(Revisit, cbind(value,100-value))
Revisit$site <- factor(Revisit$site) # just to be sure it is taken as a factor by glmmTMB
str(Revisit)

formula.MLM3.linear <- y ~ env * trait  + (1 + trait| site) + (1 + env|species)
MLM3 <- glmmTMB(formula.MLM3.linear, family = betabinomial,  data=Revisit)
summary(MLM3)

# extension with polynomial main effects

formula.MLM3.quad <- y ~ poly(env,2) + poly(trait,2) +  env : trait  + (1 + trait| site) + (1 + env|species)

MLM3.quad <- update(MLM3,formula.MLM3.quad)
summary(MLM3.quad)
anova(MLM3,MLM3.quad)
# no need for the quadratic model

# LRT test ----------
# perhaps useful if there are less than 20 species or sites
formulaNULL <- update(formula(MLM3),  ~ . - trait:env )
MLM0 <- glmmTMB(formulaNULL, data=Revisit, family= family(MLM3))
anova(MLM0,MLM3)

# Graphs of fixed and random effects ------------------------------------------------

library(ggplot2)
library(dplyr)

plot_MLM3(MLM3) +
  ggtitle('site and species effects (fixed + random) against TMG and C:N') +
  xlab("TMG                              C:N ratio") +
  ylab("estimated effects")

# LRT test ----------
# perhaps useful if there are less than 20 species or sites
formulaNULL <- update(formula.MLM3.linear,  ~ . - trait:env )
MLM0 <- glmmTMB(formulaNULL, data=Revisit, family= family(MLM3))
anova(MLM0,MLM3)


# Bootstrap confidence interval -------------------------------------------

getInteractionParameter <- function(MLM3){
  if (class(MLM3)[1] == "glmmTMB"){
    B <- summary(MLM3)$coefficients$'cond'
  } else if (class(MLM3)[1] == "lme4") {
    B <- summary(MLM3)$coefficients
  }
  coef_te <- B[nrow(B),-4]
  return(coef_te)
}
(stats_data <- getInteractionParameter(MLM3))
# not implemented in glmmTMB
# bootMer(MLM3, FUN = getInteractionParameter, nsim = 4)
#

dat_boot <- model.frame(MLM3)
nsimul <- 10 # 1000
bootstrap_mat<- matrix(0, nrow = nsimul, ncol = length(stats_data))
colnames(bootstrap_mat) <- names(stats_data)
for (i in 1:nsimul){
  dat_boot$y <- as.matrix(simulate(MLM3))
  suppressWarnings(MLM3boot <- glmmTMB(formula(MLM3), family = betabinomial,  data=dat_boot))
  bootstrap_mat[i, ] <- getInteractionParameter(MLM3boot)
}
summary(bootstrap_mat)
# 95% naive/quantile/percentile  bootstrap confidence interval (do not use..)
b_te_boot = bootstrap_mat[,1]
plot(density(b_te_boot))
b_te_qu_ci = quantile(b_te_boot,probs = c(0.025,0.975)) # 95 % percentrile bootstrap c.i.
b_te_qu_ci
mean(b_te_boot)
sd(b_te_boot)

# 95% studentized bootstap confidence interval (use this one)
b_te_est <- stats_data[1]
se_b_te_est <- stats_data[2]
se_b_te_boot <- bootstrap_mat[,2]
# t-value when the true interaction coefficient = b_te_est
tval_te_boot = (b_te_boot - b_te_est)/se_b_te_boot
plot(density(tval_te_boot))
z_b_te_qu = quantile(tval_te_boot,probs = c(0.025,0.975), na.rm=TRUE)
z_b_te_qu # ideally symmetic -2 to 2, if not, the interval corrects for bias

b_te_Stu_boot_ci = b_te_est - rev(z_b_te_qu)*se_b_te_est
names(b_te_Stu_boot_ci)=rev(names(b_te_Stu_boot_ci))
b_te_Stu_boot_ci
# Studentized intervals via boot.ci of library boot
t0 <- c(b_te_est, se_b_te_est^2)
t <- bootstrap_mat
t[,2] <- t[,2]^2
objboot <- list(t0 = t0, t = t, R = nrow(t))
boot::boot.ci(objboot, type = "stud")

# Bootstrap test ----------------------------------------------------------
set.seed(123)
nrepet <- 19
system.time(res_boot <- Bootstrap_test_prmtrc(MLM3, test_stat = "Wald", nrepet = nrepet, Binomial_total = 100))
names(res_boot)
round(res_boot$p_values,3)


# model-based permutation max test ----------------------------------------

set.seed(323)
system.time(res_perm <- MLM3_p_max(MLM3, test_stat = "Wald", nrepet = nrepet, Binomial_total = 100, print = 1))
names(res_perm)
round(res_perm$p_values,3)

