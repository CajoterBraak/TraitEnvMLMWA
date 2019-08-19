#' @title Weigthed averaging methods for pairwise trait-environment associations
#'
#' @description
#' \code{CWMSNC_regressions} carries out weighted or unweighted CWM/SNC regression for each combition of trait and environmental variable
#' with permutation testing via the max test and computation of the (un)weighted fourth-corner correlations.
#' The output can be summarized and plotted using \code{summary} and \code{plot}.
#' @param  E a data frame, numeric n-vector, factor or n-by-p matrix, containing the values of the environmentalvariables in n sites, or
#' a TE_obj object containing the trait, environment and abundance data,
#' @param  L a n by m matrix or data frame of abundance values of m species in n sites
#' @param  T a data frame, numeric m-vector, factor or m-by-q matrix, containing the values of the trait(s) for m species
#' @param  weighing type of weights used in the regression of CWM on E and SNC on T; "unw", "FC" or  "N2" (default)
#' @param  cutoff minimal number of occurrences of species
#' @param  nrepet number of permutations
#' @param  how_to_permute a list for two \code{\link[permute]{how}} calls from the permute package, the first for site permutation, the second for species permutation.
#' Default: completely random permuation.
#' @details   \code{CWMSNC_regressions} is for pairwise association between traits and environmental variables.
#' Factors in E and T data frames are expanded in to dummy variable format. For p-values of an entire factor, use \code{\link{WA_p_max}}.
#' Entering E and T as matrices or data frames has the advantage the variable names are retained.
#'
#' The permutation test(s) are carried
#' out by the max test. This test takes the maximum of the p-values of the CWM-E regression(s) and
#' the SNC-T regression(s). Both p-values are obtained by permutation.
#' In the former, the rows of E are permuted (i.e. the sites),
#' in the latter the rows of T are permuted (i.e. the species).
#'
#' Details of the  \code{weighing} options:
#' "FC" give the results of the unembellished fourth-corner correlation (which uses weighting by site and species totals);
#' "unw" gives each site and species equal weight in the regression;
#' "N2" weighs each site and species in proportion to its Hill number of order 2 (N2);
#' i.e. the effective richness of the site and effective number of occurrences of the species.
#' All weighing variants use the same CWM and SNC values (the weights in calculating CWM and SNC being the abundances in L).
#'
#' The option \code{nperm} in \code{\link[permute]{how}} is overruled by \code{nrepet}.
#' @return  An object of class CWMSNCr which is a named list, among which,
#' \item{p_values}{three permutational p-values: site-based and species-based and the maximum of these two values}
#' \item{p_parametric}{parametric p-values: site-based and species-based and the maximum of these two values
#' (can be trusted only for weighing =='unw'). NA for more than one trait or one environmental variable.}
#' \item{wFC}{(weighted) fourth-corner correlations for sites and species and its signed minimum.}
#' \item{lm_CWMe}{lm-object of the regression of CWM on E}
#' \item{lm_SNCt}{lm-object of the regression of SNC on T}
#' @references
#' ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#'
#' ter Braak, C.J.F., Peres-Neto, P.R. & Dray, S. (2018) Simple parametric tests for trait–environment association.
#' Journal of Vegetation Science, 29, 801-811. https://doi.org//10.1111/jvs.12666
#'
#' Peres-Neto, P.R., Dray, S. & ter Braak, C.J.F. (2017) Linking trait variation to the environment:
#' critical issues with community-weighted mean correlation resolved by the fourth-corner approach. Ecography, 40, 806-816.
#' https://doi.org/10.1111/ecog.02302
#' @seealso \code{\link{plot.CWMSNCr}} and \code{\link{summary.CWMSNCr}}.
#' @example demo/TutorialWA.r
#' @import permute
#' @export CWMSNC_regressions
#' @export
#'
#'
# Appendix to ter Braak 2019
#New robust weighted averaging- and model-based methods for assessing trait-environment relationships.

# functions for the statistical analysis of trait-environment association
# 1 CWMSNC_regressions: Weighted averaging (WA) regressions (CWM/SNC regressions) with permutation testing,
#                       with and without N2- or RK (Fourth-corner correlations) weighting
# 2 summary.CWMSNCr   : summary of CWMSNC_regressions
# 3 plot.CWMSNCr      : Plots of CWMSNC_regressions (main effects and trait-env association)
#
# the remaining functions are for internal use in the above three main functions
CWMSNC_regressions<-function(E, L, T, weighing = "N2", cutoff = 0, nrepet = 499, how_to_permute=list(sites=how(),species=how())){
  # combination of CWM- and SNC-based regressions, possibly weighted.
  # @ weighing = 0 or "unw": no weights
  #          = 1 or "RK" or "FC"  weights are the row and columns totals of L, as in the original fourth-corner
  #          = 2 or "N2" weights are the Hill N2-effective numbers
  #               (effective number of occurrences of a species and N2-diversity of a site)
  # @cutoff  if E is not a TE_obj, cutoff on minimal number of occurrences of a species
  dummy_coding <- function(x, dat){model.matrix(as.formula(paste("~ 0", paste(x,  collapse= "+"), sep= "+")), data = dat)}
  #end local function
  # E and T might be data frames with factors: make sure each factor has all levels included in the design matrix
  if (is.factor(E) ||is.factor(T)) {
    if (is.factor(E)) E <- as.data.frame(E)
    if (is.factor(T)) T <- as.data.frame(T)
  }
  if (is.data.frame(E)) E <- Reduce(cbind,lapply(names(E), function(x){dummy_coding(x,E)}))
  if (is.data.frame(T)) T <- Reduce(cbind,lapply(names(T), function(x){dummy_coding(x,T)}))
  obj <- make_obj_for_traitenv(E,  L,  T, cutoff = cutoff)
  result <- CWMSNC_regressions0(obj, weighing = weighing)
  # correlations
  result$wcorCWMSNC <- with(result, array(c(wcor(E, CWM),wcor(SNC,T))),dim = c(ncol(obj$E),ncol(obj$T),2),
                            dimnames = list(colnames(obj$E),colnames(obj$T), c("CWM_E_cor","SNC_T_cor")) ) # this wcor gvies E by T matrix
  result$wFC <- with(result,FC_cor_generalized(obj, wE = wsites, wT = wspecies, R = R, K = K, cutoff = cutoff, wNam = weighing))
  #names(result$wFC) <- paste(weighing, c("_sites","_species", "_fourth_corner"),"_cor", sep="")
  # p-values of association
  result$p_values <- WA_p_max1(obj, weighing= weighing, nrepet =nrepet, wsvd = FALSE, type = 'onebyone', how_to_permute=how_to_permute)$p_values
  attr(result$p_values, which = "nrepet") <- nrepet
  class(result) <- "CWMSNCr"
  return(result)
}




# the remaining functions are for internal use in the above three main functions

FC_cor_generalized <- function(E, L=0, T=0, wE=rep(1,nrow(L)) , wT= rep(1,ncol(L)),  wNam = "w", R = rowSums(L), K = colSums(L), cutoff = 0){
  # generalized fourth-corner correlation: wE and wT are weights for sites and species respectively.
  # if wE == R and wT == K it gives the fourth-corner correlation
  # value: three correlations, one for sites and one for species and the signed minimum (0 if different in sign)

  # local function
  f_cor_min <- function(x,y){
    # takes the signed minimum correlation if the correlations have the same sign (close to min)
    # 0  if the correlations differ in sign

    # if ((is.na(x)|is)) return(NA)

    # if (x*y < 0 ) { # the two elements in x diffent sign
    #   g <- 0
    # } else {
    #   g <-  sign(x[1]) * min(abs(x)) # min
    # }
    #id <- which(is.na(x)|is.na(y))
    g <- ifelse(x*y < 0, 0, sign(x)*pmin(abs(x),abs(y)))
    g
  }
  # end local function
  if(class(E)[1]!="TE_obj"){
    obj <- make_obj_for_traitenv(E,L, T,cutoff)
  } else obj <- E
  wE <- wE/sum(wE)
  wT <- wT/sum(wT)
  result <- with(obj, {
    Estd_w <- standardize_w(E,wE)
    Tstd_w <-  standardize_w(T,wT)
    Lw_si <- diag(wE/R)%*% L
    Lw_sp <-  L %*% diag(wT/K)
    #Fourthcorner  <- t(Estd_w)  %*% (L%*% Tstd_w)
    r_wE <- (t(Estd_w) %*% (Lw_si %*% Tstd_w))/sum(wE) # sum(wE)== sum(Lw_si)
    #r_wT <- sum((t(Lw_sp)%*% Estd_w) *  Tstd_w)/sum(wT) # sum(wT)== sum(Lw_sp)
    r_wT <- (t(Estd_w) %*% (Lw_sp %*% Tstd_w))/sum(wT) # sum(wT)== sum(Lw_sp)
    result <- array(c(r_wE, r_wT,f_cor_min(r_wE,r_wT)),dim = c(ncol(obj$E),ncol(obj$T),3),
                    dimnames = list(colnames(obj$E),colnames(obj$T),paste(wNam, c("FC_r_sites", "FC_r_species", "FC_r_min"), sep = "_")))

    return(result)
  })
  return(result)
}



# Hill number of order 2: N2
fN2 <- function(x){x <- x/sum(x); 1/sum(x*x)}


CWMSNC_regressions0<-function(E,L=0,T=0, weighing = 0, cutoff = 0){
# combination of CWM- and SNC-based regressions, possibly weighted.
# @param E  object from class TE_obj from make_obj_for_traitenv(E,L,T,cutoff) or n-vector or nx1 matrix with the environmental values
# @param L  n x m matrix of abundance values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
# @param T  m-vector or m x 1 matrix with trait values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
# @ weighing = 0 or "unw": no weights
#          = 1 or "RK" or "FC" weights are the row and columns totals of L, as in the original fourth-corner
#          = 2 or "N2" weights are the Hill N2-effective numbers
#               (effective number of occurrences of a species and N2-diversity of a site)
# @cutoff  if E is not a TE_obj, cutoff on minimal number of occurrences of a species

  if(class(E)[1]!="TE_obj"){
    obj <- make_obj_for_traitenv(E, L = L,T = T,cutoff)
  } else obj <- E
  result <- with(obj,{
      R <- rowSums(L) # the site totals
      K <- colSums(L) # the species totals
      if (weighing %in% c(0,"unw") ){
        wsites <- rep(1, nrow(L))
        wspecies <-rep(1,ncol(L))
        weighing = "unw"
      } else if (weighing %in% c(1,"RK","FC")){
        wsites <- R
        wspecies <-K
        weighing = "RKw"
      } else if (weighing %in% c(2,"N2")){
        N2_si <- apply(L, 1, fN2)
        N2_sp <- apply(L, 2, fN2)
        wsites <- N2_si
        wspecies <- N2_sp
        weighing = "N2w"
      }
      CWM <- (L%*%T)/R #Community weighted mean wrt to trait T (mean trait value per site)
      SNC <- (t(L)%*%E)/K  # Species niche centroid wrt to environmental variable E (mean environmental value per species)
      meanE <- mean_w(E)
      meanT <- mean_w(T)
      wmeanCWM <- mean_w(CWM, w = wsites)
      wmeanSNC <- mean_w(SNC, w = wspecies)
      # site-level regression of CWM on E
      lm_CWMe <- lm(CWM~E, weights = wsites)
      if (ncol(CWM) ==1) {
        anova_site <- anova(lm_CWMe)
        F_sites <- anova_site$`F value`[1]# get F-value and p-value of site-level test
        Prob.site<- anova_site$"Pr(>F)"[1]
      } else {F_sites <-NA; Prob.site <- NA}
      # species-level regression of SNC on T
      lm_SNCt <- lm(SNC~T, weights =   wspecies)
      if (ncol(SNC)==1){
        anova_species <- anova(lm_SNCt)
        F_species <- anova_species$`F value`[1]# get F-value and P-value of species-level test
        Prob.species <- anova_species$"Pr(>F)"[1]
      } else {F_species <-NA; Prob.species <- NA}
      result<-list(p_parametric = c(p.site.prmtrc=Prob.site,p.species.prmtrc=Prob.species,p.max.prmtrc=max(Prob.site,Prob.species)), F_values =  c(F_sites = F_sites, F_species= F_species),
                   lm_CWMe = lm_CWMe, lm_SNCt = lm_SNCt, CWM = CWM, SNC = SNC,
                   E = E, T = T, L = L, meanE = meanE, meanT = meanT, wmeanCWM = wmeanCWM, wmeanSNC = wmeanSNC,
                   R=R, K= K, wsites = wsites, wspecies=wspecies,weighing=weighing,cutoff=cutoff)
    return(result)
    }
  )
  return(result)
}

CWMr<-function(E,L=0,T=0, weighing = "unw", cutoff = 0){
  # CWMr CWM-based regressions, possibly weighted.
  # @param E  object from class TE_obj from make_obj_for_traitenv(E,L,T,cutoff) or n-vector or nx1 matrix with the environmental values
  # @param L  n x m matrix of abundance values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
  # @param T  m-vector or m x 1 matrix with trait values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
  # @ weighing = 0 or "unw": no weights
  #          = 1 or "RK" or "FC" or "regular" weights are the row and columns totals of L
  #          = 2 or "N2" weights are the Hill N2-effective numbers
  #               (effective number of occurrences of a species and N2-diversity of a site)
  # @cutoff  if E is not a TE_obj, cutoff on minimal number of occurrences of a species

  if(class(E)[1]!="TE_obj"){
    obj <- make_obj_for_traitenv(E,L,T,cutoff)
  } else obj <- E
  result <- with(obj,{
    R <- rowSums(L) # the site totals
    if (weighing %in% c(0,"unw") ){
      wsites <- NULL
    } else if (weighing %in% c(1,"RK", "FC")){
      wsites <- R
    } else if (weighing %in% c(2,"N2")){
      N2_si <- apply(L, 1, fN2)
      wsites <- N2_si
    }
    CWM <- (L%*%T)/R #Community weighted mean wrt to trait T (mean trait value per site)
    # site-level regression of CWM on E
    lm_CWMe <- lm(CWM~E, weights = wsites)
    anova_site <- anova(lm_CWMe)
    F_sites <- anova_site$`F value`[1]# get F-value and p-value of site-level test
    Prob.site<- anova_site$"Pr(>F)"[1]
    result<-list(p_values =Prob.site, F_value = F_sites,
                 lm_CWMe = lm_CWMe,  CWM = CWM,  E = E, T = T, L = L, wsites = wsites,weighing=weighing,cutoff=cutoff)

    return(result)
  }
  )
  return(result)
}


SNCr<-function(E,L=0,T=0, weighing = 0, cutoff = 0){
  # CWMr CWM-based regressions, possibly weighted.
  # @param E  object from class TE_obj from make_obj_for_traitenv(E,L,T,cutoff) or n-vector or nx1 matrix with the environmental values
  # @param L  n x m matrix of abundance values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
  # @param T  m-vector or m x 1 matrix with trait values or nonspecified if E = object from function make_obj_for_traitenv(E,L,T,cutoff)
  # @ weighing = 0: no weights
  #          = 1 or "RK" or "FC" or  weights are the row and columns totals of L
  #          = 2 or "N2" weights are the Hill N2-effective numbers
  #               (effective number of occurrences of a species and N2-diversity of a site)
  # @cutoff  if E is not a TE_obj, cutoff on minimal number of occurrences of a species

  if(class(E)[1]!="TE_obj"){
    obj <- make_obj_for_traitenv(E,L,T,cutoff)
  } else obj <- E
  result <- with(obj,{
    #R <- rowSums(L) # the site totals
    K <- colSums(L) # the species totals
    if (weighing %in% c(0,"unw") ){
      wsites <- NULL
      wspecies <-NULL
    } else if (weighing %in% c(1,"RK", "FC")){
      #wsites <- R
      wspecies <-K
    } else if (weighing %in% c(2,"N2")){
      N2_sp <- apply(L, 2, fN2)
      wspecies <-N2_sp
    }
    SNC <- (t(L)%*%E)/K  # Species niche centroid wrt to environmental variable E (mean environmental value per species)
    #species-level regression of SNC on T
    lm_SNCt <- lm(SNC~T, weights =   wspecies)
    anova_species <- anova(lm_SNCt)
    F_species <- anova_species$`F value`[1]# get F-value and P-value of species-level test
    Prob.species<- anova_species$"Pr(>F)"[1]
    result<-list(p_values =Prob.species, F_value = F_species,
                 lm_SNCt = lm_SNCt,  SNC = SNC,  E = E, T = T, L = L, wspecies = wspecies,weighing=weighing,cutoff=cutoff)

    return(result)
  }
  )
  return(result)
}

CWMSNCr <- function(obj, options){
  # level in c("species","sites","both")
  if(class(obj)[1]!="TE_obj") {print("obj must be of class TE_obj in CWMSNCr"); return(obj)}
  WA <- with(options,{
  if (level == "sites") WA <- c(site_p = CWMr(obj, weighing = weighing)$p_values) else
    if (level == "species") {
      WA <- c(species_p =SNCr(obj, weighing =weighing )$p_values)
    } else if (level == "both") {
      WA <- c(CWMr(obj, weighing =weighing)$p_values,SNCr(obj, weighing =weighing )$p_values)
      names(WA)<- c("site_p", "species_p")
    }
  return(WA)
  })
  if (options$invert_p_value) WA <- 1/WA
  return(WA)
}

WA_p_max0 <- function(obj, nrepet = 19,  weighing = "N2") {
  #n_sites <- nrow(obj$L); n_species <- ncol(obj$L)

  options1 <- setoptions4pmax_test(verbose = FALSE,
                       FUN_test_statistics = list(init=function(obj,options){list(obj =obj, options = options)}, test=CWMSNCr),
                       with.MLM= FALSE, nrepet = nrepet,
                       perm.mat.spe = 0 , perm.mat.sit = 0, print = 1000,
                       nmax = 6, filek = -1)
  options2 <- setoptions4WA(weighing = weighing , invert_p_value = TRUE, fast = FALSE)
  options3 <- c(options1,options2[-length(options2)])
  result <- PermutationTest_r_c_max(obj, options = options3)
  return(result)
}



# local functions
mean_w <- function(X,w = rep(1/nrow(X),nrow(X))){t(w/sum(w))%*% X}
center_w <- function(X,w = rep(1/nrow(X),nrow(X))){ X - rep(1,length(w))%*%t(w)%*% X }
standardize_w <- function(X,w = rep(1/nrow(X),nrow(X)), wsvd = FALSE){
  # NB requires w to be have sum 1
  ones <- rep(1,length(w))
  Xc <- X - ones %*% t(w)%*% X
  Xstd <- Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w))
  if (wsvd) Xstd <- wSVD(Xstd, w)
return(Xstd)
}

wcor <- function(X, Y=X, w = rep(1,nrow(X))){
  # weighted correlation between matrix X and Y
  w <- w/sum(w)
  Xstd <- standardize_w(X, w)
  Ystd <- standardize_w(Y, w)
  t(Xstd) %*% diag(w) %*% Ystd
  }

wSVD <- function(Y,w=rep(1/nrow(Y),nrow(Y))){
  sw <- sqrt(w)
  Ystar <- Y*sw
  svdY <- svd(Ystar)
  Ustar <- svdY$u
  id <- which(svdY$d > 1.e-6)
  return(Ustar[,id, drop = FALSE]/sw)
} # returns w-orthogonalized Y

init_WA4_pmax <- function(obj, options){
 res <- CWMSNC_regressions0(obj, weighing = options$weighing)
 Wn <- res$wsites/sum(res$wsites)
 Ws <- res$wspecies/sum(res$wspecies)
 E <-standardize_w(res$E, Wn, wsvd = options$wsvd)
 T <- standardize_w(res$T, Ws, wsvd = options$wsvd)
 CWM <- center_w(res$L%*%T/res$R, Wn) #CWM wrt to standardized T (trait), and with wSVD orthogonalized, centered for numerical precision
 SNC <- center_w(t(res$L)%*%E/res$K, Ws) # SNC wrt to standardized E (environment) and with wSVD orthogonalized, centered for numerical precision
 result <- list(E=E,T=T, CWM=CWM, SNC = SNC, p_parametric = res$p_parametric, Wn = Wn, Ws=Ws)
 return(list(obj = result, options = options))
}

#obj_init <- init_WA4pmax(obj,options)

test_stat_WA <- function(obj, options){
  test_stat <- with(obj, {
    if (options$level =="sites"){
      Estd <- standardize_w(E, Wn, wsvd = options$wsvd)
      test_stat <- t(Estd)%*%diag(Wn)%*% CWM # E by T
      #test_stat <- sum(test_stat*test_stat)
      if (options$type == 'joint') {
        test_stat2 <-t(rep(1,ncol(Estd))) %*% (test_stat*test_stat)%*% rep(1,ncol(CWM))
      } else {
        test_stat2 <- c(abs(test_stat))
      }
    } else {
      Tstd <- standardize_w(T, Ws, wsvd = options$wsvd)
      test_stat <- t(SNC)%*% diag(Ws)%*% Tstd  # t(Tstd)%*% diag(Ws)%*% SNC
      if (options$type == 'joint') {
       test_stat2 <- t(rep(1,ncol(SNC))) %*% (test_stat*test_stat)%*% rep(1,ncol(Tstd))
      } else {
        test_stat2 <- c(abs(test_stat))
      }
      return(test_stat2)
    }
     #return(sum(test_stat*test_stat))
    })
 return(test_stat)
}

#test_stat_WA(obj_init,options)
WA_p_max1 <- function(obj, nrepet = 19,  weighing = "N2", type= 'joint', wsvd = TRUE, how_to_permute=list(sites=how(nperm = nrepet),species=how(nperm = nrepet))) {
  #n_sites <- nrow(obj$L); n_species <- ncol(obj$L)

  options1 <- setoptions4pmax_test(verbose = FALSE,
                                   FUN_test_statistics = list(init= init_WA4_pmax, test= test_stat_WA),
                                   with.MLM= FALSE, nrepet = nrepet,
                                   perm.mat.spe = 0 , perm.mat.sit = 0, print = 5000,
                                   nmax = 6, how_to_permute=how_to_permute, filek = -1)
  options2 <- setoptions4WA(weighing = weighing, invert_p_value = FALSE, type = type, wsvd = wsvd, fast = TRUE)
  options3 <- c(options1,options2[-length(options2)])
  if (options3$type != 'joint') options3$wsvd <- FALSE
  result <- PermutationTest_r_c_max(obj, options = options3)
  # reform result in case of type separate
  if (options3$type != "joint"){ #onebyone
     result$p_values <- array(result$p_values, dim = c(ncol(obj$E),ncol(obj$T),3),
                  dimnames = list(colnames(obj$E),colnames(obj$T), c("p.sites","p.species","p.max")))
  }
  return(result)
}


PermutationTest_r_c_max <- function(obj, options = c(setoptions4pmax_test(),setoptions4MLM()),  ...) {
  # Specialized max test for WA-based and model-based MLM based methods
  #function to determine randomization p-values based on test statistics for rows and columns
  # based on row and column based permutations of (possibly model-based) test-statistics,
  # such as the anova  F-value or chisq-value or equivalently, 1/p-values
  # The function examimes exceedance and takes the absolute value of the test-statistics, that explains why 1/p-value can be used.
  #
  # @param obj object of class TE_obj (from make_obj_for_traitenv)
  # @param in options FUN_test_statistics  function that returns one or more test statistics (value, vector or matrix with one of two rows)
  #                        it should return 1/pval, because the
  #
  # e.g.FUN_test_statistic =
  # ... options for FUN_test_statistic
  # BEWARE: the size calculation is adapted from Miller et al 2019 MEE https://doi.org//10.1111/2041-210X.13119 but
  # is only trustworthy if the null model holds true!!!!
  # For a non-null model the size is biased upwards (ter Braak, 2019 unpubl.)
  #
  # For  permutation of model-based methods trait0 and env0 are added to obj, these are not permuted, whereas E and T are.
  # the model formula are therefore adapted to contain the terms env0 and trait0


  control_perm <- options$how_to_permute
  if (is.matrix(options$perm.mat.spe)) options$perm.mat.spe <- options$perm.mat.spe[sample(nrow(options$perm.mat.spe)),]
  if (is.matrix(options$perm.mat.site)) options$perm.mat.site <- options$perm.mat.site[sample(nrow(options$perm.mat.site)),]

  result <- options$FUN_test_statistics$init(obj,options)
  obj <- result$obj; options <- result$options
  if (!is.null(obj$p_parametric))p_parametric = obj$p_parametric else p_parametric <- NA
  rm(result)
  options_sites <- options_species <- options
  options_sites$verbose  <- options_species$verbose <- FALSE
  options_species$level <- "species"
  options_sites$level <- "sites"
  if (options$with.MLM){
    # add random T and E to data: trait0 and env0 and update formulas accordingly
    i_spe0 <- sample(1:nrow(obj$T), 1)
    i_sit0 <- sample(1:nrow(obj$E), 1)
    per.col <- syssample(i = i_spe0, options$perm.mat.spe, n =ncol(obj$L),control_perm = control_perm[[2]])
    obj$trait0 <- matrix(obj$T[per.col,,drop =FALSE],nrow= nrow(obj$T), ncol =ncol(obj$T))  # used in species level
    per.row <- syssample(i = i_sit0, options$perm.mat.site, n =nrow(obj$L), control_perm = control_perm[[1]])
    obj$env0 <- matrix(obj$E[per.row,,drop =FALSE],nrow= nrow(obj$E), ncol =ncol(obj$E)) # used in site level

    #  e.g:  y ~  trait*env  + (1 + trait| site) + (1 + env|species) becomes
    #       y ~  trait*env  + trait0 + (1 + trait| site) + (1 + env|species) #for species
    # #     y ~  trait*env  + env0 + (1 + trait| site) + (1 + env|species) # for sites

    options_species$formula <-  update(options_species$formula, ~ . + trait0 )
    options_sites$formula  <-  update(options_sites$formula, ~ . + env0 )
  }
  obs.species <- options_species$FUN_test_statistics$test(obj, options = options_species)
  obs.sites <- options_sites$FUN_test_statistics$test(obj, options = options_sites)
  obs <- rbind(obs.sites,obs.species)
  nrepet <- options$nrepet
  if (options$verbose) {cat("test statistics for the data\n"); print(obs); cat("progress out of ", nrepet, " permutations:\n")}
  if (options$with.MLM){
    # change trait0 and env0 to T and E in data and update formulas so that only env:trait is permuted
    obj$trait0 <- obj$T
    obj$env0 <- obj$E
    #  e.g:  y ~  trait*env  + trait0 + (1 + trait| site) + (1 + env|species) becomes
    #       y ~  trait*env  + trait0 + (1 + trait0| site) + (1 + env|species) #for species
    # #     y ~  trait*env  + env0 + (1 + trait| site) + (1 + env0|species) # for sites

    options_species$formula <- update(options_species$formula,~ . -(1+trait|site)+(1+trait0|site))
    options_sites$formula <- update(options_sites$formula,~ . -(1+env|species)+(1+env0|species))

  }
  if(is.matrix(options$perm.mat.spe)|| is.matrix(options$perm.mat.sit)){
    nrepet <- min(c(nrow(options$perm.mat.spe)-1,nrow(options$perm.mat.sit)-1, nrepet) )#  minus identity ispe_0 isit_0 permutation (which gives collinearity)
  } else { nrepet <- options$nrepet}
  sim.row <- matrix(0, nrow = nrepet, ncol = ncol(obs))
  sim.col <- matrix(0, nrow = nrepet, ncol = ncol(obs))


  for(i in 1:nrepet){
    per.row <- syssample(i = i, exclude = i_sit0, perm.mat = options$perm.mat.site, n =nrow(obj$E), control_perm = control_perm[[1]])
    per.col <- syssample(i = i, exclude = i_spe0, perm.mat = options$perm.mat.species, n =nrow(obj$T), control_perm = control_perm[[2]])
    #permute_rows_columns <- function(obj, per.row,per.col){
    obj.row <- obj.col <-  obj  # obj.rc <- obj
    obj.row$E <-obj.row$E[per.row,,drop =FALSE]
    obj.col$T <-obj.col$T[per.col,,drop =FALSE]
    suppressWarnings(sim.row[i, ] <- options_sites$FUN_test_statistics$test(obj.row,options = options_sites))
    suppressWarnings(sim.col[i, ] <- options_species$FUN_test_statistics$test(obj.col,options = options_species))
    if (!i%%options$print){
      iteration <- i
      if ((i/options$print)%%10 && i != nrepet) cat(iteration, " ") else cat(iteration, "\n ")
      if (options$filek >=0 ) save(options, sim.row, sim.col, obs, iteration, file =paste("perm_r_c",options$filek,".rdata", sep = ""))
    }
  }
  ialpha <- 1/options$alpha
  if (ncol(obs)==1){
    isna.r <-  sum(is.na(sim.row))
    isna.c <-  sum(is.na(sim.col))
    pval.row <- (sum(abs(sim.row) >= abs(obs[1,]), na.rm = TRUE) + 1)  / (nrepet- isna.r  + 1)
    pval.col <- (sum(abs(sim.col) >= abs(obs[2,]), na.rm = TRUE) + 1)  / (nrepet- isna.c  + 1)
    if (options$falseSize){
      size.row <- (sum(sim.row >= ialpha, na.rm=TRUE))/ (nrepet - isna.r)
      size.col <- (sum(sim.col >= ialpha, na.rm=TRUE))/ (nrepet - isna.c)
      size.max <- (sum(max(sim.row,sim.col) >= ialpha, na.rm=TRUE))/ (nrepet - isna.c)
    }
  } else {
    obs.mat.row <- matrix(rep(abs(obs[1,]),each=nrepet), nrow= nrepet, ncol=ncol(obs))
    isna.r <-  colSums(is.na(sim.row))
    pval.row <- (colSums(abs(sim.row) >= obs.mat.row, na.rm=TRUE) + 1)/ (nrepet - isna.r + 1)
    isna.c <-  colSums(is.na(sim.col))
    obs.mat.col <- matrix(rep(abs(obs[2,]),each=nrepet), nrow= nrepet, ncol=ncol(obs))
    pval.col <- (colSums(abs(sim.col) >= obs.mat.col, na.rm=TRUE) + 1)/ (nrepet - isna.c + 1)
    if(options$falseSize){
      size.row <- (colSums(sim.row >= ialpha, na.rm=TRUE))/ (nrepet - isna.r)
      size.col <- (colSums(sim.col >= ialpha, na.rm=TRUE))/ (nrepet - isna.c)
      colnames(sim.row) = colnames(obs)
      colnames(sim.col) = colnames(obs)
    }
  }
  result <- c(p.site.permut = pval.row, p.species.permut = pval.col, pmax.permut = pmax(pval.row, pval.col))
  attr(result, "nrepet")<- nrepet
  if (options$falseSize) result <- cbind(result, t(rbind(size.row = size.row, size.col = size.col)))
  return(list(p_values=result, nrepet = nrepet, obs = obs,sim.row=sim.row, sim.col = sim.col, p_parametric = p_parametric))
}


syssample <- function(i = 0, exclude = 0, perm.mat = 0, n, control_perm = permute::how()){
  # function to generate either random samples (if !is.matrix(perm.mat)) or
  # a systematic sample: i-th or for i>= exclude: i+1 th row of permat (i.e. without the excluded sample)
  # exclude = 0 for no exclude, or permutation number to exclude
  if (is.matrix(perm.mat)){
    if (!exclude){ ii <- i} else{
      if(i<exclude) ii <- i else {ii <- i + 1}
    }
    sam <- perm.mat[ii,]
  } else {
    sam <-permute::permute(i= i, n=n, control = control_perm)
  }
  return(sam)
}

f_allPerms <- function(n){
  # all permutation matrix without identity and reverse order
  result <- permute::allPerms(n)[1:(numPerms(n)/2 - 1)  ,]
  result[sample(nrow(result)),]
}

