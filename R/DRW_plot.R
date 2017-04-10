# DRW package - plot state

#' Plot DRW state on go
#'
#' @param mpdt
#' mobile particle data table (x, y, L)
#' @param rpdt
#' source release data table (x, y, L)
#' @param drts
#' integer [1]; DRW time step
#' @param mfsp
#' integer [1]; MODFLOW stress period
#' @param bas
#' BAS.MFpackage
#' @param wel
#' WEL.MFpackage or NULL
#' @param gccs,grcs
#' numeric []; column and row divider co-ordinates
#' @param ...
#' graphical parameters such as \code{xlim} (not axis labels, titles,
#'  \code{asp} or \code{zlim})
#'
#' @return
#' \code{NULL}
#'
#' @import data.table
#' @import graphics
#' @importFrom Rflow MFimage
#'
plotDRWstate <- function(mpdt, rpdt, drts, mfsp,
                         bas, wel, gccs, grcs, ...){
  for(plL in unique(mpdt$L)){
    # plot MODFLOW model domain
    MFimage(bas$IBOUND[,, plL], gccs, grcs, c(-1, 1),
            c("blue", "grey", "transparent"), asp = 1, ...,
            xlab = "x", ylab = "y",
            main = paste0("layer ", plL, "\ntimestep ", drts),
            sub = paste(if(is.null(mpdt)) "0" else nrow(mpdt),
                        "mobile particles"))

    # plot mobile particle swarm
    if(!is.null(mpdt) && nrow(mpdt) > 1L)
      points(mpdt[L == plL & m != 0, .(x, y)],
             pch = 16L, cex = .4, col = "#00000080")

    # plot active source releases
    if(!is.null(rpdt) && nrow(rpdt) > 1L)
      points(rpdt[L == plL & mtmp != 0, .(x, y)],
             pch = 16L, cex = .8, col = "purple")

    # plot wells if available
    if(!is.null(wel)){
      # active wells in this layer
      points(wel$data[L == plL & sp == mfsp & Q != 0, {
        list(x = mean(gccs[C + 0:1]),
             y = mean(rev(grcs)[R + 1:0]))
      }, by = c("C", "R")][, .(x, y)],
      pch = 10L, cex = 2, lwd = 2, col = "red")

      # active wells in other layers
      points(wel$data[L != plL & sp == mfsp & Q != 0, {
        list(x = mean(gccs[C + 0:1]),
             y = mean(rev(grcs)[R + 1:0]))
      }, by = c("C", "R")][, .(x, y)],
      pch = 10L, cex = 1, lwd = 1.5, col = "darkred")
    }
  }

  invisible(NULL)
}
