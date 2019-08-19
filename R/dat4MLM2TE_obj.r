#' @title Utility function dat4MLM2TE_obj
#'
#' @description
#' \code{dat4MLM2TE_obj} creates from data in a format for glm(m)/MLM a \code{TE_obj} object (see \code{\link{make_obj_for_traitenv}}), which consists of
#' trait, environment and abundance data. The argument \code{dat} is usually created using
#' \code{\link{expand4glmm}}. The function is used repeatedly in \code{\link{MLM3_p_max}}.
#' @param  dat dataframe with names y, site, species, trait and env
#' for values of abundance, identity of site and of species and the trait value and environmental value, respectively.
#' @inheritParams make_obj_for_traitenv
#' @details
#' The order of data in dat should be in standard order, that is, first all data from site 1, then site 2,
#' with all species listed in identical order in each site (including zero abundances).
#' BEWARE: The function currently works for a single trait and environmental variable only!!!!
#' @return  An object of class TE_obj
#' @examples
#' \dontrun{
#' data("Revisit")
#' str(Revisit)
#' TE_obj <- dat4MLM2TE_obj(Revisit)
#' str(TE_obj)
#' }
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @export dat4MLM2TE_obj
dat4MLM2TE_obj<- function(dat, cutoff = 0){
  # create abundance matrix TE_obj from a dataframe dat that is in standard order
  #     first all data from site 1, then site 2,
  #        with the all species listed in each site (so include 0 abundances) and
  #        species having the same order in each site
  # @param dat dataframe with names y, sites, species
  # create abundance matrix L, environmental variable vector E and trait vector T
  # BEWARE: currently works for a single trait and environmental variable only!!!!
  dat$site=factor(dat$site)
  dat$species=factor(dat$species)
  n_sites <-nlevels(factor(dat$site))
  n_species <-nlevels(dat$species)
  species <- dat$species[seq(from  = 1, by = n_sites, length.out = n_species)]
  sites <- dat$site[1:n_sites]
  # create the usual trait T, environment E and abundance matrix L;
  # these would normally the original data.. but are here derived from the existing data.frame dat
  T <- matrix(dat$trait[seq(from  = 1, by = n_sites, length.out = n_species)], nrow = n_species, dimnames = list(species, "trait"))
  E <- matrix(dat$env[1:n_sites], nrow = n_sites, dimnames = list(sites, "env"))
  if (length(dat$y)>nrow(dat)){ # e.g. with binomial response
    L <- matrix(dat$y[,1], nrow = n_sites,ncol = n_species, dimnames = list(sites= sites,species=species))
  }
  else { # e.g. with count data
    L <- matrix(dat$y, nrow = n_sites,ncol = n_species, dimnames = list(sites= sites,species=species))
  }
  # check the data for empty species and sites and bring them in standard format of class TE_obj
  obj <- make_obj_for_traitenv(E = E, L = L,T = T, cutoff = cutoff)
}
