#' calculate reversed cumulative sum
#' @param a vector
#' @return a vector of reversed cumulative sum
#' @examples cumsum_rev(1:10)
#' 55 54 52 49 45 40 34 27 19 10
#'
#'
cumsum_rev<- function(a) rev(cumsum(rev(a)))

#' calculate sum of rows
#' @param x a matrix to calculate rowsums
#' @return a vector of row sums
#'
#'
rowSums <- function (x, na.rm = FALSE, dims = 1L)
{
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2L)
    return(x)
  if (dims < 1L || dims > length(dn) - 1L)
    stop("invalid 'dims'")
  p <- prod(dn[-(id <- seq_len(dims))])
  dn <- dn[id]
  z <- if (is.complex(x))
    .Internal(rowSums(Re(x), prod(dn), p, na.rm)) + (0+1i) *
    .Internal(rowSums(Im(x), prod(dn), p, na.rm))
  else .Internal(rowSums(x, prod(dn), p, na.rm))
  if (length(dn) > 1L) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[id]
  }
  else names(z) <- dimnames(x)[[1L]]
  z
}

#' transform Hasse diagram to paritial order inequalities
#' @param Hasse the partial order relations expressed by directed graph
#' @param Z a matrix of staging variable (poset)
#' @return a matrix expressing the patial order constraints inequalities on coefficients
#' @import igraph

as_inequalities <- function(Hasse, Z){
  cellName = colnames(Z)
  PO = matrix(0, nrow = ecount(Hasse), ncol = length(cellName))
  colnames(PO) = cellName
  edgelist = as_edgelist(Hasse)
  for (i in 1:ecount(Hasse)){
    PO[i,edgelist[i,1]] = -1
    PO[i,edgelist[i,2]] = 1
  }
  return(PO)
}


