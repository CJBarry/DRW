# DRW package - sorption and desorption

#' Sorb and Desorb
#'
#' @param mpdt,ipdt data.table (ts, x, y, L, zo, m);
#' for the mobile and immobile states respectively
#' @param Rf numeric;
#' retardation factor, conceptually intended to be >= 1
#'
#' @return
#' named list of data.tables for the new mobile and immobile states
#'
#' @import data.table
#'
sorb.desorb <- function(mpdt, ipdt, Rf){
  newi <- mpdt[, list(ts = ts, x = x, y = y,
                       L = L, zo = zo, m = m*(1 - 1/Rf))]
  mpdt[, m := m/Rf]

  # sorbed mass giving back
  newm <- ipdt[, list(ts = ts, x = x, y = y,
                       L = L, zo = zo, m = m/Rf)]
  ipdt[, m := m*(1 - 1/Rf)]

  list(mob = rbind(mpdt, newm),
       immob = rbind(ipdt, mewi))
}