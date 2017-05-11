# DRW package - kernel smoothed concentration

#' Weighted Kernel Density Estimate of Concentration
#'
#' A post-processing tool for \code{\link{DRW}}.  This function gives a
#'  distributed concentration estimate for a DRW result, using a weighted
#'  kernel density estimate from the particles, using
#'  \code{\link[ks]{kde}}.
#'
#' @param DRWmodel
#' \code{DRWmodel} S4 object or a character string giving the file path to
#'  one
#' @param mfdata,wtop
#' open NetCDF, character string file path(s) or list of open NetCDFs;
#' MODFLOW data set(s) and saruated water tops
#' @param dkcell
#' numeric [1] or [2];
#' cell spacing for the kernel smooth in x and y directions
#' @param ts,L
#' \code{NULL} or integer[];
#' time steps/ layers to process (\code{NULL} implies all)
#' @param smd
#' numeric [1];
#' smoothing distance for kernel smooth: too small and the result will be
#'  patchy, too large and the result will be too spread
#' @param H
#' matrix [2, 2];
#' alternative, lower level way of specifying smoothing (see
#'  \code{\link[ks]{kde}})
#' @param binned
#' logical [1];
#' see \code{\link[ks]{kde}}; \code{NA} uses binning for particle sets
#'  numbering more than 500 and not otherwise
#' @param Kxlim
#' \code{"auto"}, \code{"mf"}, numeric [2] or something that, with
#'  \code{Kylim}, can be read by \code{\link[grDevices]{xy.coords}};
#' \code{Kxlim} and \code{Kylim} give the x and y limits for the region in
#'  which to determine the concentration; "auto" uses the particle
#'  distribution and "mf" uses the MODFLOW model edges
#' @param Kylim
#' \code{NULL} or numeric [2]
#' @param ...
#' additional arguments for \code{\link[ks]{kde}}
#' @param ptype
#' \code{"plume"} (the default) or \code{"sorbed"}
#'
#' @return
#' List object of class kDRW with elements:\cr
#' \code{$conc} numeric [x, y, L, ts]; concentration estimate\cr
#' \code{$coords} list:\cr
#' \code{..$x,y} numeric []; x and y grid divide co-ordinates (note not the
#'  cell centres)\cr
#' \code{..$L,ts} integer []; layers and time steps that are present\cr
#' \code{$H} numeric [2, 2]; the H matrix used for smoothing (see
#'  \code{\link[ks]{kde}})
#'
#' @importFrom ks kde
#' @import Rflow
#' @import RNetCDF
#' @import data.table
#' @import methods
#' @importFrom grDevices xy.coords
#' @importFrom stringr str_dup
#' @export
#'
#' @examples
#' setwd(mfdir <- system.file(package = "DRW"))
#' drmod <- readRDS("DRW_EXAMPLE.rds")
#'
#' kdr <- ksConc(drmod, "drw_mf_demo.nc", "drw_mf_demo_wtop.nc", 10,
#'               smd = 20, Kxlim = "mf")
#'
#' # Note here that the result is indexed using strings for the 3rd and 4th
#' #  dimensions.  This is much less likely to cause confusion because the
#' #  layers and time steps included are likely not to be a simple sequence
#' #  starting from 1.
#' with(kdr, image(z = conc[,, "L1", "ts10"], coords$x, coords$y,
#'                 col = grDevices::rainbow(51L)))
#'
ksConc <- function(DRWmodel, mfdata, wtop, dkcell, ts = NULL, L = NULL,
                   smd, H, binned = NA,
                   Kxlim = c("auto", "mf")[1L], Kylim = NULL,
                   ..., ptype = "plume"){
  # get time values (if given) and particle data table
  DRWmodel <- switch(class(DRWmodel)[1L],
                     character = readRDS(DRWmodel),
                     DRWmodel = DRWmodel,
                     stop("ksConc: invalid DRWmodel"))
  tvals <- DRWmodel@time
  pdt <- slot(DRWmodel, ptype)

  # get mfdata list
  mfdatal <- switch(class(mfdata)[1L],
                    character = {
                      l <- lapply(mfdata, open.nc)
                      on.exit(lapply(mfdatal, close.nc), add = TRUE)
                      l
                    },
                    NetCDF = list(mfdata),
                    list = mfdata,
                    stop("ksConc: invalid mfdata"))
  wtopl <- switch(class(wtop)[1L],
                  character = {
                    l <- lapply(wtop, open.nc)
                    on.exit(lapply(wtopl, close.nc), add = TRUE)
                    l
                  },
                  NetCDF = list(wtop),
                  list = wtop,
                  stop("ksConc: invalid wtop"))
  #
  # - get data set timing information
  #  -- time steps
  mftstl <- lapply(mfdatal, mftstime, absolute = TRUE)

  #
  #  -- MODFLOW model start times and final end time
  dsett <- c(lapply(mftstl, `[`, 1L), last(last(mftstl)), recursive = TRUE)

  # check consistency between mfdata and DRWmodel if possible
  if(is(DRWmodel, "DRWmodel")){
    if(!identical({
      MFinfo <- sapply(c("title", "author", "history"),
                       function(att){
                         sapply(mfdatal, att.get.nc,
                                "NC_GLOBAL", att)
                       })
      if(is.vector(MFinfo)) MFinfo <- t(MFinfo)
      MFinfo
    }, DRWmodel@MFinfo))
      warning("ksConc: mfdata does not match the MODFLOW data sets used for DRWmodel")
  }

  # determine which rows of pdt to use
  rtu <- (if(is.null(ts)) TRUE else pdt$ts %in% ts) &
    (if(is.null(L)) TRUE else pdt$L %in% L)
  pdt <- pdt[rtu]
  Ls <- sort(unique(pdt$L))
  tss <- sort(unique(pdt$ts))

  # get x and y limits
  if(is.character(Kxlim)){
    Kxlim <- switch(Kxlim[1L],
                    auto = list(x = unname(quantile(pdt$x,
                                                    c(.01, .99), TRUE)),
                                y = unname(quantile(pdt$y,
                                                    c(.01, .99), TRUE))),
                    mf = list(x = range(gccs(mfdatal[[1L]], TRUE)),
                              y = range(grcs(mfdatal[[1L]], TRUE))),
                    stop("ksConc: invalid Kxlim"))
  }
  xylims <- xy.coords(Kxlim, Kylim)
  Kxlim <- xylims$x; Kylim <- xylims$y
  #
  # - adjust x and y limits so that their ranges are exact multiples of
  #    dkcell values
  dkcell <- switch(length(dkcell), rep(dkcell, 2L), dkcell)
  if(is.null(dkcell))
    stop("ksConc: invalid dkcell (should be length 1 or 2 numeric)")
  #
  #  -- determine mismatches
  xmm <- diff(Kxlim)%%dkcell[1L]
  ymm <- diff(Kylim)%%dkcell[2L]
  #
  #  -- adjust symmetrically
  if(xmm > 0) Kxlim <- Kxlim + (dkcell[1L] - xmm)*c(-.5, .5)
  if(ymm > 0) Kylim <- Kylim + (dkcell[2L] - ymm)*c(-.5, .5)

  # set up results array
  # - x and y grid divider co-ordinates
  xgcs <- seq(Kxlim[1L], Kxlim[2L], dkcell[1L])
  ygcs <- seq(Kylim[1L], Kylim[2L], dkcell[2L])
  #
  # - eval.points
  xpts <- (xgcs[-1L] + xgcs[-length(xgcs)])/2
  ypts <- (ygcs[-1L] + ygcs[-length(ygcs)])/2
  epts <- expand.grid(x = xpts, y = ypts)
  setDT(epts)
  #
  #  -- determine C and R references and water thickness at eval points
  epts[, C := cellref.loc(x, gccs(mfdatal[[1L]], TRUE))]
  epts[, R := cellref.loc(y, grcs(mfdatal[[1L]], TRUE), TRUE)]
  thk <- array(0, c(length(xpts), length(ypts), length(Ls), length(tss)),
               list(NULL, NULL, paste0("L", Ls), paste0("ts", tss)))
  epts[, {
    for(Ln in Ls){
      bot <- nc.imtx(mfdatal[[1L]], "elev",
                     cbind(C, R, Ln + 1L))

      for(tsn in tss){
        ds <- pdt[ts == tsn, unique(mfds)]
        mfts <- pdt[ts == tsn, unique(mfts)]
        stopifnot(length(ds) == 1L)
        stopifnot(length(mfts) == 1L)

        top <- nc.imtx(wtopl[[ds]], "wtop", cbind(C, R, Ln, mfts))
        HDRY <- att.get.nc(mfdatal[[ds]], "Head", "HDRY")
        top[top == HDRY] <- NA_real_

        thk[,, paste0("L", Ln), paste0("ts", tsn)] <<- top - bot
      }
    }
    NULL
  }]
  #
  # - the array
  #  -- the third and fourth dimensions have names to refer to the correct
  #      layers and timesteps (L3, L4, L5, ... and ts6, ts7, ts8, ... for
  #      example)
  ar <- array(0, c(length(xpts), length(ypts), length(Ls), length(tss)),
              list(NULL, NULL, paste0("L", Ls), paste0("ts", tss)))

  # fill results array with results
  if(missing(H)) H = diag(smd^2L, 2L)
  autobinned <- is.na(binned)[1L]
  cat(str_dup(" ", nch <- 3L + 4L + 5L + 4L + 4L))
  pdt[, {
    k <- kde(cbind(x, y), H = H, eval.points = epts[, list(x, y)],
             binned = if(autobinned) .N > 500L else binned,
             w = m/mean(m), ...)

    # k$estimate gives values for normalised mass per 2D volume (or area)
    # to get mass in 2D cell:
    mkcell <- k$estimate*prod(dkcell)*sum(m)

    # volume in each k cell
    Vkcell <- prod(dkcell)*thk[,, paste0("L", L), paste0("ts", ts)]

    ar[,, paste0("L", L), paste0("ts", ts)] <<- mkcell/Vkcell

    cat(str_dup("\b", nch), "L: ", formatC(L, width = 4L),
        str_dup(" ", 5L), "ts: ", formatC(ts, width = 4L), sep = "")
    NULL
  }, by = c("L", "ts")]; cat("\n")

  # return result and meta data
  structure(list(conc = ar,
                 coords = list(x = xgcs, y = ygcs, L = Ls, ts = tss),
                 H = H),
            class = "kDRW")
}
