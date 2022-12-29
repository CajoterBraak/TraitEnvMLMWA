#' @title Plotting the interaction of a MLM3 model
#'
#' @description
#' \code{plot_MLM3} plots the trait-environment interaction
#' of a fitted MLM3 object.
#' @param  MLM3 the fitted MLM3 object, created by glmer (lme4) or glmmTMB.
#' @param  verbose logical; default FALSE. If TRUE, a summary of the MLM3 object is printed.
#' @param  title character text to adapt the main title.
#' @details  The code uses libraries ggplot2 and dplyr. The trait and environment variable names
#' in the MLM3 model must be 'trait' and 'env' as in the example (as they are hard coded)
#' @return  A list of printable and modifiable objects of class ggplot2.
#' @seealso \code{\link{expand4glmm}}.
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @examples
#' \dontrun{
#' data("Revisit")
#' formula.MLM3 <- y ~  env*trait  + (1 + env|species) + (1 + trait| site)
#' MLM3 <- glmmTMB(formula.MLM3, family = betabinomial,  data=Revisit)
#' plot_MLM3(MLM3)
#' }
#' @import ggplot2
#' @export
plot_MLM3 <- function(MLM3, verbose = FALSE, title = paste( 'site and species effects (fixed + random) against env and trait' )){
  # plot of a fitted MLM3 object
  # the names in the data.frame (model.frame) of MLM3 should be in any order
  # y, species, sites, trait, env
  # data should be in column-wise order (first all the data of the first species, then the second....)
  # get data in Y, trait, env form (Figure 1)----------------------------------------------------------------
  dat <-model.frame(MLM3)
  n_sites <- with(dat, nlevels(factor(site)))
  n_species <- with(dat, nlevels(factor(species)))
  species <- dat$species[seq(from  = 1, by = n_sites, length.out = n_species)]
  sites <- dat$site[1:n_sites]
  trait <- dat$trait[seq(from  = 1, by = n_sites, length.out = n_species)]
  env <- dat$env[1:n_sites]
  if (length(dat$y) > nrow(dat)){
    # binomial data
     Y <- matrix(dat$y[,1], nrow = n_sites,ncol = n_species, dimnames = list(sites= sites,species=species))
  } else {
    # count data
     Y <- matrix(dat$y, nrow = n_sites,ncol = n_species, dimnames = list(sites= sites,species=species))

  }

  # extract random effect and b_te non-null ------------------------------------------

  if (verbose) print(summary(MLM3))
  if (class(MLM3)[1] == "glmmTMB"){
    B <- summary(MLM3)$coefficients$'cond'
  } else { # lme4
    B <- summary(MLM3)$coefficients
  }
  b_te <- B[nrow(B),1]; gamma0 <- B["env",1]; beta0 <- B["trait",1]
  xx <- ranef(MLM3)$'cond'
  b_i <- xx$site$trait; c_j <-xx$species$env

  b_i_star <- beta0 + b_te*env + b_i      # Eq (2a)
  c_j_star <- gamma0 + b_te*trait + c_j   # Eq (2b)

  ETdat <- data.frame(env, b_i_star, N = rowSums(Y>0), level = "sites")
  names(ETdat)[2]<- "b_star"
  ETdat$beta0 = beta0
  ETdat$b_tslope = beta0 + b_te*env

  TEdat <- data.frame(trait, c_j_star= c_j_star, N = colSums(Y>0), level = "species")
  names(TEdat)[2]<-"c_star"
  TEdat$gamma0 = gamma0
  TEdat$b_eslope = gamma0 + b_te*trait

  ETTEdat <- suppressMessages(dplyr::full_join(ETdat, TEdat, by= c("N","level")))
  # facet plot --------------------------------------------------------------
  pp <- ggplot(ETTEdat) +
    geom_point(aes(x= env, y = b_star, size = N), shape = 1) +
    geom_point(aes(x= trait, y = c_star, size = N), shape = 16) +
    geom_line(aes(x = env, y = b_tslope), size = 2) +
    geom_line(aes(x = trait, y = b_eslope), size = 2) +
    geom_hline(yintercept = 0, linetype = 'dotted',size = 1) +
    geom_line(aes(x = env, y = beta0), linetype = 'longdash', size = 1) +
    geom_line(aes(x = trait, y = gamma0), linetype = 'longdash', size = 1) +
    facet_wrap(~level, scales = "free") +
    ggtitle(title) +
    xlab("environmental variable                                    trait") +
    ylab("estimated effects")
  return(pp)
}
