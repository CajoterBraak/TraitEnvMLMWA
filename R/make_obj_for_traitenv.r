#' @title Utility function make_obj_for_traitenv
#'
#' @description
#' \code{make_obj_for_traitenv} combines the three data items needed for analyzing trait-environment association
#' into a single object of class \code{TE_obj}. In the process, empty sites and non-occurring species are deleted. There is an option
#' for stronger deleting criteria.
#' @param  E a numeric n-vector or n by p matrix containing p environmental variable(s) of n sites
#' @param  L a n by m matrix or data frame of abundance values of m species in n sites
#' @param  T a numeric m-vector or m by q matrix containing q traits of m species
#' @param  cutoff minimal number of occurrences of species
#' @return  A named list,
#' \item{L}{a sites by species matrix of abundances, after removal of species and sites, based on cutoff}
#' \item{E}{a matrix of environmental values (rows are the same sites as the rows of L) }
#' \item{T}{a matrix of trait values (rows are the same as the column of L)}
#' @details  If E or T is a matrix, the column names are retained. If E or T is a vector, E gets name "env" and T gets name "trait".
#' If you wish to retain the original name,
#' enter E and T as data frames or single column matrices with a column name.
#' @references ter Braak (2019) New robust weighted averaging- and model-based methods
#' for assessing trait-environment relationships. Methods in Ecology and Evolution (https://doi.org/10.1111/2041-210X.13278)
#' @seealso \code{\link{expand4glmm}} and \code{\link{WA_p_max}}.
#' @examples
#' \dontrun{
#' data("aravo",  package = "ade4")
#' Y <-aravo$spe
#' trait <- scale(aravo$traits$SLA)
#' env <- scale(aravo$env$Snow)
#' obj <- make_obj_for_traitenv(env,Y, trait, cutoff=0)
#' str(obj)
#' }
#' @export make_obj_for_traitenv

make_obj_for_traitenv <- function(E, L, T, cutoff = 0){
  # makes object of class TE_obj from matrices or dataframes from L,E and T (L=Y)
  # no factors allowed, as all is converted to matrices
  # the object is a list of three matrices: L, E, T
  # @param cutoff species occurring strictly less then cutoff are deleted

  L <- as.matrix(L)
  E <- as.matrix(E)
  T <- as.matrix(T)
  if (nrow(E)!= nrow(L)){
    stop("number of rows in Env E not equal to number of rows of species abundance matrix L")
  } else if (nrow(T)!= ncol(L)){
    stop("number of rows in Trait T not equal to number of columns of species abundance matrix L")
  }
  rows<-seq_len(nrow(L))
  cols<-seq_len(ncol(L))
  rni <-which(rowSums(L)==0)
  repeat {
    if (length(rni)) {L <- L[-rni,,drop = FALSE]; rows <-rows[-rni]}
    ksi <- which(colSums(L)==0)
    if (length(ksi)) {L <- L[,-ksi, drop = FALSE]; cols <- cols[-ksi]}
    rni <-which(rowSums(L)==0)
    if ( length(rni)==0 & length(ksi)==0){break}
  }
  E <-as.matrix(E)[rows,,drop = FALSE]
  T <-as.matrix(T)[cols,,drop = FALSE]
  # end check_L()
  if(cutoff >0 ){
    abOcc <- apply(L>0,2,mean)
    L <- L[,abOcc>cutoff, drop =FALSE]
    T <- as.matrix(T)[abOcc>cutoff,, drop = FALSE]
  }
  #print(str(L))
  if(is.null(rownames(L))) rownames(L)<- paste("site", 1:nrow(L), sep="")
  if(is.null(colnames(L))) colnames(L)<- paste("spec",1:ncol(L),sep="")
  rownames(E)= rownames(L); if(is.null(colnames(E))){ if (ncol(E)>1) colnames(E) <-  paste("E", 1:ncol(E), sep = "") else colnames(E) <- "env"}
  rownames(T)= colnames(L); if(is.null(colnames(T))){ if (ncol(T)>1) colnames(T) <-  paste("T", 1:ncol(T), sep = "") else colnames(T) <- "trait"}
  obj = list(L=L, E=E,T = T)
  class(obj) = c("TE_obj", class(obj))
  return(obj)
}
