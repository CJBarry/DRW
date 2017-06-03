# DRW package - source term

#' DRW source term from a DNAPL source term
#'
#' see \code{\link[DNAPL]{DNAPLSourceTerm}}
#'
#' @param dnst
#' \link[DNAPL]{DNAPLSourceTerm} object
#' @param mfdata
#' list of NetCDFs of MODFLOW output
#' @param gccs,grcs
#' numeric[];
#' column and row divider co-ordinates
#'
#' @import RNetCDF
#' @import data.table
#' @importFrom stats punif
#' @importFrom stats approxfun
#'
#' @return
#' data.table, source term for \code{\link{DRW}}
#'
ST.DNAPL <- function(dnst, mfdata, gccs, grcs){
  hL <- dnst@DNAPLmodel@hL
  tmp <- data.table(x = dnst@xy[1L],
                    y = dnst@xy[2L],
                    z = dnst@z0 + rev(cumsum(rev(hL))) - hL/2)
  tmp[, C := cellref.loc(x, gccs)]
  tmp[, R := cellref.loc(y, grcs, TRUE)]

  # correct for NA z0, assuming that the base of the model is desired
  if(is.na(dnst@z0)){
    tmp[, z := {
      base <- var.get.nc(mfdata[[1L]], "elev",
                         c(C, R,
                           dim.inq.nc(mfdata[[1L]], "NLAY")$length + 1L),
                         c(1L, 1L, 1L))
      base + rev(cumsum(rev(hL))) - hL/2
    }, by = c("C", "R")]
  }

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
    apply(dnst@Jeffluent, 1L, approxfun,
          x = dnst@time, yleft = 0)
  }]
  tmp[, c("C", "R", "z") := NULL]
  tmp
}
