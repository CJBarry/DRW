# DRW package - coalescing

#' Dynamic Coalescing
#'
#' @param pdt
#' particle data.table (ts, x, y, L, zo, m)
#' @param cd
#' numeric [2];
#' coalescing search radii in the horizontal and vertical directions
#' @param mm
#' numeric [1];
#' minimum mass of particles after coalescing
#' @param maxnp
#' integer [1] or \code{Inf};
#' maximum number of particles after coalescing
#' @param mfdata
#' NetCDF object; MODFLOW data set
#' @param wtop
#' NetCDF object;
#' cell by cell top of water (see \code{\link[Rflow]{get.wtop.nc}})
#' @param mfts
#' integer [1];
#' current MODFLOW time step number
#'
#' @return
#' list with:\cr
#' \code{$state}: particle data.table with column and row references
#'  (ts, x, y, L, zo, m)\cr
#' \code{$loss}: mass lost due to being in inactive cells
#'
#' @import coalesce
#' @import data.table
#' @import Rflow
#' @importFrom stats weighted.mean
#'
coalesceDRW <- function(pdt, cd, mm, maxnp, mfdata, wtop, mfts){
  # in the case that all mass is depleted, return a zero-row data table
  if(all(pdt$m == 0)) return(data.table(ts = integer(0L),
                                        x = double(0L),
                                        y = double(0L),
                                        L = integer(0L),
                                        zo = double(0L),
                                        m = double(0L)))

  # The treatment the vertical dimension in this coalescing algorithm is a
  #  bit different.  Rather than using the z co-ordinate, the z-offset
  #  within layers is used, in order to be more true to the way MODPATH
  #  works, and to the geological structure.  Firstly, coalescence is done
  #  layer by layer, so that no layer can lose mass to another layer during
  #  this process.  And within layers, the centre of mass will be preserved
  #  in terms of z-offset, but not necessarily in terms of z.  z will not
  #  be perfectly preserved in cases when mass moves across a boundary
  #  between two cells for which the saturated top and bottom are
  #  different.  This also means that the loss of mass from the top and
  #  bottom of the model is impossible, although loss into the corner of an
  #  inactive cell is still possible, though rare.  The vertical search
  #  radius should still be given in units of length.  This function will
  #  find an appropriate conversion to a z-offset-based search radius based
  #  on the weighted average saturated thickness of the cells occupied by
  #  particles in each layer

  # find column and row references for each particle
  # - probably, in fact, these have already been found from the lost mass
  #    step before coalescing
  if(!"C" %in% names(pdt)){
    pdt[, C := cellref.loc(x, gccs(mfdata, TRUE))]
  }
  if(!"R" %in% names(pdt)){
    pdt[, R := cellref.loc(y, grcs(mfdata, TRUE), TRUE)]
  }

  # find layer thicknesses
  pdt[, thk := if(is.na(C) || is.na(R) || is.na(L)) NA else{
    top.imtx <- cbind(C, R, L, mfts)
    bot.imtx <- cbind(C, R, L + 1L)

    top <- nc.imtx(wtop, "wtop", top.imtx)
    bot <- nc.imtx(mfdata, "elev", bot.imtx)

    top - bot
  }, by = c("C", "R", "L")]

  # register particle mass that is in inactive cells and discard from state
  #  data table
  loss <- sum(pdt[is.na(thk), m])
  pdt <- pdt[!is.na(thk)]
  if(!nrow(pdt)) return(list(state = pdt[, .(ts, x, y, L, zo, m)],
                             loss = loss))

  # coalesce layer by layer
  copdt <- pdt[, {
    # weighted average thickness
    wat <- weighted.mean(thk, m)

    # determine appropriate z-offset-based vertical search radius based on
    #  cd[2L] and wat
    cdzo <- cd[2L]/wat

    # perform coalesce
    tmp <- coalesce(.SD[, .(x, y, zo, m)], cd[1L], cdzo, mm, maxnp = maxnp,
                    subregions = TRUE, nppsr = 256L, TwoD = FALSE)
    tmp[, ts := ts[1L]]
    tmp
  }, by = L]
  setcolorder(copdt, c("ts", "x", "y", "L", "zo", "m"))

  list(state = copdt, loss = loss)
}
