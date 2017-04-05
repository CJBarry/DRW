# DRW package - utility functions

#' Index matrix for NetCDF
#'
#' Emulate the \code{A[B]} matrix indexing behaviour in R for NetCDF
#'  connections whilst avoiding reading large arrays into memory.
#'
#' @param ncfile
#' NetCDF object (see \code{\link[RNetCDF]{open.nc}})
#' @param variable
#' character string; the variable name in the NetCDF file
#' @param imtx
#' integer matrix;
#' the indices to read; must have the same number of columns as variable
#'  has dimensions; must contain only positive integers and no NAs; ideally
#'  the range of values in some or most columns would not be very large to
#'  avoid reading a large array into working memory
#'
#' @return
#' a vector whose type corresponds to that of \code{variable}
#'
#' @import RNetCDF
#'
nc.imtx <- function(ncfile, variable, imtx){
  stopifnot(is(ncfile, "NetCDF"))
  stopifnot(is.character(variable))
  stopifnot(is.matrix(imtx))

  ndims <- as.integer(var.inq.nc(ncfile, variable)$ndims)

  # check that the index has the correct number of columns
  stopifnot(identical(ndims, ncol(imtx)))

  # find index ranges for each dimension
  rgs <- apply(imtx, 2L, range)

  # read a sub-array from the NetCDF, with only the required dimension
  #  ranges
  # - this is a compromise between speed and memory: the smallest possible
  #    array that would require only one read from the NetCDF
  ar <- var.get.nc(ncfile, variable, rgs[1L,], apply(rgs, 2L, diff) + 1L,
                   collapse = FALSE)

  imtx_mod <- vapply(1:ndims, function(d) imtx[, d] - rgs[1L, d] + 1L,
                     integer(nrow(imtx)))
  if(is.vector(imtx_mod)) imtx_mod <- t(imtx_mod)

  ar[imtx_mod]
}
