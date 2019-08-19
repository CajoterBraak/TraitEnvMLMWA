#' @title Summary of a weigthed averaging-based analysis of trait-environment association
#'
#' @description
#' \code{summary.CWMSNCr} summarizes the results of \code{\link{CWMSNC_regressions}}. For a single quantitative trait
#' and a single quantitative environmental variable three types of correlations and p-values are
#' given (sites, species, min/max). For nominal and multipe
#' trait and environment data, one set of such correlations and p-values is printed
#' (by default, the max p-values and signed minimum fourth-corner correlations)
#' with an associated heatmap of the significant associations. The p-values can be adjusted for
#' multiple testing.
#' @param  object an object of class CWMSNCr, created by \code{\link{CWMSNC_regressions}}.
#' @param  digits number of digits to print
#' @param  type  "sites", "species"  or "max" (default, which combines sites and species). Used only for multiple or nominal trait and environmental data
#' @param p_value_adjust_method method for adjustment of the p-values for multiple comparison.
#' "fdr" (default), "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", or "none". Used only for multiple or nominal trait and environmental data
#' @param significance_level threshold for plotting a heatmap of the signed p-values
#' @param ... other optional arguments
#' @return  For quantitative single trait-single environmental variable data:
#' a numeric matrix with trait-environment correlations (first row) and p-values (second row) for sites, species and their min/max combination (columns).
#' Correlations in the min/max column are signed minima of sites and species, being 0 if the correlations for sites and species differ in sign.
#' p-values in the min/max column are the maximum of the p-values for sites and species.
#'
#' For nominal or multiple trait/environment data: a list of
#' \item{p_val_adj}{the adjusted p-values}
#' \item{FCcorrelations}{the (weighted) fourth corner correlations}
#' \item{heatmap}{a gglot object containing the basis heatmap on the basis of argument \code{significance_level}}
#' @details \code{p_value_adjust_method} uses  \code{\link[stats]{p.adjust}} on P-values of all
#' pxq trait-environment combinations of the site-level and the species-level tests,
#' which, for \code{type == 'max'} are subsequently combined by taking the maximum of
#' the adjusted site P-value of the CWM regression and adjusted species P-value for the SNC regression
#' for each combination.
#' This adjustment method is the same as used in \code{\link[ade4]{fourthcorner}} and is less
#' conservative than applying p.adjust on the max P-values directly.
#' @seealso \code{\link{plot.CWMSNCr}}.
#' @example demo/TutorialWA_Aravo.r
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @export
summary.CWMSNCr <- function(object, ..., digits = 3, type ="max", p_value_adjust_method = "fdr", significance_level = 0.05){
  # object should be an object resulting from CWMSNC_regressions
  if (sum(dim(object$wFC)[-3])==2){ # single trait and env analysed
    res <- rbind(object$wFC, object$p_values)
    rownames(res) <- paste(c("correlations", "p-values"), " (", object$weighing, ")" , sep = "")
    colnames(res) <- c("sites","species","min/max")
    attr(res, "nrepet")<- attr(object$p_values, which = "nrepet")
    print(round(res, digits = digits))
    invisible(res)
  } else { # more than one trait or environmental variable
    itype <- switch (type,
      sites = 1,
      species = 2,
      max = 3
    )
      if (itype ==3) typeFC = "signed min" else typeFC <- type
      p_val_adj <- matrix(p_adjust(object$p_values, p_value_adjust_method, itype), nrow= dim(object$p_values)[1], ncol =dim(object$p_values)[2], dimnames = dimnames(object$p_values)[c(1,2)])
      cat(paste("Fourth-corner correlations", " (", object$weighing, ")", " (", typeFC, ")\n" , sep = ""))
      print(round(object$wFC[,,itype], digits = digits))
      cat(paste("p-values", " (", object$weighing, ")", " (", type, ")\n" , sep = ""))
      print(round(p_val_adj, digits = digits))
      sp_values <- p_val_adj * ifelse(sign(object$wFC[,,itype])==0,1,sign(object$wFC[,,itype]))
      # print p-values only with sign of the correlation
      #  print(round(sp_values, digits = digits))
      hm <- heatmap_p(sp_values, significance_level = significance_level)
      print(hm)
      result <- list(p_val_adj =p_val_adj, FCcorrelations = object$wFC[,,itype], heatmap = hm)
      invisible(result)
  }
}

p_adjust <- function(p_values, p_value_adjust_method, itype){
  #p_values: array of ExTx c("sites", "species" and "max") for itype = 1,2,3
  if (itype == 3){
    # Adjust site and species separately, as above and then use maxP
    p.sites <- p.adjust(p_values[,,1],method = p_value_adjust_method)
    p.species <- p.adjust(p_values[,,2],method = p_value_adjust_method)
    p_val_adj <- pmax(p.sites,p.species)
  } else if (itype == 1){
   # Adjust p-values of CWMr regressions
    p_val_adj <- p.adjust(p_values[,,1],method = p_value_adjust_method)
  } else { #if (itype == 2){
    # Adjust p-values of SNCr regressions
    p_val_adj <- p.adjust(p_values[,,2],method = p_value_adjust_method)
  }
  return(p_val_adj)
}

heatmap_p <- function(x, significance_level = 0.05){
  # heatmap using ggplot
  xx<-x
  xx[abs(x) > significance_level] <- 0
  xx[ xx < 0] <- -1
  xx[ xx > 0] <- 1
  df <- expand.grid(Env = rownames(x), Traits=colnames(x))
  df$association = factor(c(xx), levels=c(-1,0,1), labels = c("-","0","+"))
  hp <- ggplot(data = df, aes(x= Env, y = Traits)) + geom_tile(aes(fill = association)) +
    theme(axis.text.x=element_text(angle=90,hjust=1)) #+ scale_fill_brewer("RdYlGn")
  hp
}


