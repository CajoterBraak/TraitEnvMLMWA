#' @title Plotting the results of a weigthed averaging-based analysis of pairwise trait-environment association(s)
#'
#' @description
#' \code{plot.CWMSNCr} plots the results of \code{\link{CWMSNC_regressions}} for the single quantitative trait
#' and the single quantitative environmental variable,
#' specified by the arguments \code{trait} and \code{env}.
#' Both main effects (site totals vs E and species totals vs T) and the associations (CWM vs E and SNC vs T) are plotted.
#' @param  x an object of class CWMSNCr, created by \code{\link{CWMSNC_regressions}}.
#' @param  title  character text for the main title of the association graph (optional)
#' @param ... other optional arguments
#' @param trait trait name or level of a nominal tait  to be plotted.
#' @param env environmental variable name or level of a nominal environmental variable to be plotted
#' @details
#' In the one E vector - one T vector case, you can choose your own labels for \code{trait} and \code{env}.
#' In the  E or T matrix or data frame case, names shoud refer to names of variables or labels of factors.
#' If in error, a list of names is returned.
#' The arguments \code{trait} and \code{env} can be numbers referring to the
#' columns numbers in x$E and x$T;
#' these numbers do not reflect the orginal E and T if any variable in them is a factor (nominal).
#'
#' The code uses the libraries ggplot2 and dplyr.
#' @return  A list of two printable and modifiable objects of class ggplot2. The first plots the main effects, the second the association.
#' The y-coordinate of the points in the first plot is the logarithm of the site total (left) and the species total (right).
#' The y-coordinate of the points in the second plot is weighted trait mean (CWM, left) and weighted environmental mean (SNC, right).
#' The long-dash line is at the weigthed mean of the y-coordinate of the points in plots. The dotted line is at the unweighted mean of the trait (left)
#' and on of the environmental variable (right). If, in the second-left plot, the dashed line is above the dotted line, the trait main effect is positive
#' (check in first plot, right diagram). If, in the second-right plot, the dashed line is above the dotted line, the environmental main effect is positive
#' (check in first plot, left diagram). In this sense, the second plot shows the sign of the main effects, shown more explicitly in the first plot.
#' @seealso \code{\link{CWMSNC_regressions}}.
#' @example demo/TutorialWA_Aravo.r
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @import ggplot2
#' @export
plot.CWMSNCr <- function(x, ..., title = NA, trait = 1, env = 1){
    # plot of a CWMSNCr object
  if (sum(dim(x$wFC)[-3])==2){ # special case one trait one env
    if (is.character(env)) {namE <- env} else {namE <- colnames(x$E)}
    if (is.character(trait)) {namT <- trait; tait <- 1} else {namT <- colnames(x$T); tait <- 1}
    env <- i <- 1
    trait <- j <- 1
  } else {
     if(is.character(env) | is.character(trait)  & (!env %in% colnames(x$E) | !trait %in% colnames(x$T))){
       if (is.character(trait) & !trait %in% colnames(x$T)) {
         message("variable ", trait, " not in trait data. Possible names are: ", paste(colnames(x$T), collapse = " "))
         message("Possible env names are: ", paste(colnames(x$E), collapse = " "))
         return(list(env.names = colnames(x$E), trait.names =colnames(x$T)))
       }
       if (is.character(env) & !env %in% colnames(x$E)) {
          message("variable ", env, " not in environmental data. Possible names are: ", paste(colnames(x$E), collapse = " "))
          message("Possible trait names are: ", paste(colnames(x$T), collapse = " "))
          return(list(env.names = colnames(x$E), trait.names =colnames(x$T)))
        }
    }
    if (is.character(env)){
      namE <- env
      i <- which(colnames(x$E)  %in% env )
    } else if (is.numeric(env)){
       namE <- colnames(x$E)[env]
       i <- env
    }
    if (is.character(trait)){
      namT <- trait
      j <- which(colnames(x$T) %in%  trait)
    } else if (is.numeric(trait)){
      namT <- colnames(x$T)[trait]
      j <- trait
    }
  }
  ETdat <- with(x, data.frame( E[,i],  CWM[,j], R, wsites, meanT[,j], wmeanCWM[,j], level = "sites"))
  names(ETdat) <- c("E", "CWM", "R", "N2", "meanT","wmeanCWM" ,"level")

  TEdat <- with(x, data.frame( T[,j],  SNC[,i], K, wspecies, meanE[,i], wmeanSNC[,i], level = "species"))
  names(TEdat) <- c("T", "SNC", "K", "N2", "meanE","wmeanSNC","level")


  ETTEdat <- suppressMessages(dplyr::full_join(ETdat,TEdat))
  if (is.na(title)[1]) {title <-paste(x$weighing, '-regressions of CWM on ', namE, " and SNC on ", namT ,collapse = "",sep="")}
  # facet plot --------------------------------------------------------------
  suppressMessages(p.main <-  ggplot(data= ETTEdat) +
                     geom_point(aes(x= E, y = log(R), size = N2),shape = 1) +  stat_smooth(aes(x= E, y = log(R), weight=N2), method = "gam") +
                     geom_point(aes(x= T, y = log(K), size = N2), shape = 16) + stat_smooth(aes(x= T, y = log(K), weight=N2), method = "gam") +
                     xlab(paste(namE, "                                                     ", namT, collapse = ""))+
                     ylab(" log totals") +
                     facet_wrap(~level, scales = "free") +
                     ggtitle(paste("log site and species totals vs", namE, "and" , namT, collapse = " ")))
  suppressMessages(p.assoc <-  ggplot(data= ETTEdat) +
                     geom_point(aes(x= E, y = CWM, size = N2),shape = 1) +  stat_smooth(aes(x= E, y = CWM, weight=N2), method = "gam") +
                     geom_point(aes(x= T, y = SNC, size = N2), shape = 16) + stat_smooth(aes(x= T, y = SNC, weight=N2), method = "gam") +
                    # xlab(" env                                                     trait")+
                     xlab(paste(namE, "                                                     ", namT, collapse = ""))+
                     ylab("weighted means") +
                     geom_line(aes(x = E, y = meanT), linetype = 'dotted', size = 1)+
                     geom_line(aes(x = T, y = meanE), linetype = 'dotted', size = 1) +
                     geom_line(aes(x = E, y = wmeanCWM), linetype = 'longdash', size = 1)+
                     geom_line(aes(x = T, y = wmeanSNC), linetype = 'longdash', size = 1) +
                     facet_wrap(~level, scales = "free") +
                     ggtitle(title))

  #  print(p.main)
  #  print(p.assoc)
  return(list(p.main,p.assoc))
}
