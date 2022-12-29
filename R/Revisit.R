#' @title Revisit data
#'
#'
#'
#' @description  The data is from the revisit to the low-elevation, non-serpentine sites of
#' Robert Whittaker’s historic plant community study sites in the Siskiyou Mountains of Southwest Oregon
#' (Damschen, Harrison & Grace 2010).
#'  Briefly, the data consists of the abundance of 75 species in 52 sites,
#'  calculated from the number of 100 quadrat corners per site that each species intersected
#'  (Miller et al. 2018, 2019).
#'  The environmental variable and functional trait are standardized versions of
#'  the original topographic moisture gradient (TMG) and the leaf carbon-to-nitrogen ratio (C:N),
#'  respectively (Miller et al. 2019). The dataframe \code{Revisit} has been derived from the
#'  whittaker_revisit_data.csv file from Miller et al. (2018), who describe the variables as follows
#' \itemize{
#' \item\code{site}  study site ID (integer, one study plot per site)
#' \item\code{species} name of species
#' \item\code{trait} functional trait: leaf tissue carbon to nitrogen ratio
#' \item\code{env} environmental variable: site position on the topographic moisture gradient (higher numbers are topographically drier sites)
#' \item\code{value} plant abundance (see Description for details)
#' \item\code{y} matrix with two columns with, as successes and failures, value and 100 - value, respectively.
#' }
#' @references
#' Damschen, E.I., Harrison, S. & Grace, J.B. (2010) Climate change effects on an endemic-rich edaphic flora: resurveying Robert H. Whittaker's Siskiyou sites (Oregon, USA). Ecology, 91, 3609-3619.
#' https://doi.org/10.1890/09-1057.1
#'
#' Miller JED, Damschen EI, Ives AR (2018) Data from: Functional traits and community composition: a comparison among community-weighted means, weighted correlations, and multilevel models. Dryad Digital Repository. https://doi.org/10.5061/dryad.7gj0s3b
#'
#' Miller JED, Damschen EI, Ives AR (2019) Functional traits and community composition: a comparison among community-weighted means, weighted correlations, and multilevel models. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210x.13119
#'
#' ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)

#' @name Revisit
#' @docType data
NULL
