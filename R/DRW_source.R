# DRW package - source term

#' DRW source term from a DNAPL source term
#'
#' see \code{\link[DNAPL]{DNAPLSourceTerm}}
#'
#' @param dnst
#'
#' @import RNetCDF
#' @import data.table
#' @importFrom stats punif
#'
#' @return
#' data.table, source term for \code{\link{DRW}}
#'
ST.DNAPL <- function(dnst){
  hL <- dnst@DNAPLmodel@hL
  tmp <- data.table(x = dnst@xy[1L],
                    y = dnst@xy[2L],
                    z = dnst@z0 +
                      rev(cumsum(rev(hL))) - hL/2)
  tmp[, C := cellref.loc(x, gccs)]
  tmp[, R := cellref.loc(y, grcs, TRUE)]
  tmp[, L := {
    cellref.loc(z,
                rev(var.get.nc(mfdata[[1L]], "elev",
                               c(C, R, NA),
                               c(1L, 1L, NA))),
                TRUE)
  }, by = c("C", "R")]
  tmp[, zo := {
    punif(z,
          var.get.nc(mfdata[[1L]], "elev",
                     c(C, R, L + 1L), c(1L, 1L, 1L)),
          var.get.nc(mfdata[[1L]], "elev",
                     c(C, R, L), c(1L, 1L, 1L)))
  }, by = c("C", "R", "L")]
  tmp[, J := {
    apply(dnst@Jeffluent, approxfun,
          x = dnst@time, yleft = 0)
  }]
  tmp[, c("C", "R", "z") := NULL]
  tmp
}
