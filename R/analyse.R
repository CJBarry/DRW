# post-processing analysis for DRW outputs

#' @rdname wellanalyse
#' @name DRWwell-analysis
#' @title Analyse receptors in a DRW model output
#'
#' @description
#' Extract the flux to individual receptors and calculate the abstracted
#'  concentration to receptors in the output of a \code{\link{DRW}} model.
#'
#' @param dr
#' a \link{DRWmodel} S4 object,
#' @param welref
#' data.table (or object coercible to one) with named columns C, R and L
#'  giving the CRL references of receptors to be analysed, integer [3],
#'  \code{NULL}, WEL.MFpackage (see \code{\link[Rflow]{read.WEL}}), or (for
#'  \code{wellJ} only) various other MFpackage objects;
#' \code{NULL} implies all cell references to which a flux has been recorded
#'  by \code{\link{DRW}}
#' @param WEL
#' WEL.MFpackage or NetCDF object of MODFLOW outputs;
#' if a NetCDF object, there must be a "Wells" array in the data set
#' @param DIS
#' DIS.MFpackage (see \code{\link[Rflow]{read.DIS}});
#' only required if \code{WEL} is a WEL.MFpackage object
#' @param MFt0
#' numeric [1];
#' only required if \code{WEL} is a WEL.MFpackage object
#' start time of the MODFLOW model, which should be the same as for the
#'  NetCDF MODFLOW data set used for the \code{\link{DRW}} model
#'
#' @return
#' \code{wellJ}: list of fluxes to each receptor specified in welref, with
#'  a value for each time step in \code{dr} (\code{dr@time})\cr
#' \code{wellC}: as for \code{wellJ} but concentrations, determined by
#'  diluting the flux with the transient abstraction rate
#'
#' In each case, the elements of the list are named giving the cell
#'  reference as a character string (e.g. \code{"C30R20L1"})
#'
NULL

#' @rdname wellanalyse
#' @import data.table
#' @export
#'
wellJ <- function(dr, welref){
  welref <- extract.welref(welref, dr, "wellJ")

  # some times the C, R and L columns are duplicated - when I sort out this
  #  minor bug, the following line can be removed:
  dr@fluxout <- dr@fluxout[, list(ts, C, R, L, J_out)]

  J <- rep(list(double(length(dr@time))), nrow(welref))
  names(J) <- welref[, paste0("C", C, "R", R, "L", L)]

  # join data.tables
  fluxes <- dr@fluxout[welref, on = c("C", "R", "L")]

  fluxes[!is.na(ts), {
    J[[paste0("C", C, "R", R, "L", L)]][ts] <<- sum(J_out)
  }, by = c("ts", "C", "R", "L")]

  # J should always be in the same order as welref
  J
}

#' @rdname wellanalyse
#' @import data.table
#' @import Rflow
#' @import RNetCDF
#' @export
#'
wellC <- function(dr, WEL, DIS, MFt0, welref = WEL){
  # check consistency of arguments according to required class combinations
  if(is(WEL, "WEL.MFpackage")){
    if(!is(DIS, "DIS.MFpackage"))
      stop("wellC: DIS package needed for stress period timings; see Rflow::read.DIS")

    if(missing(MFt0))
      stop("wellC: MFt0 required to give global start time of MODFLOW model")

    if(is(MFt0, "NetCDF"))
      MFt0 <- att.get.nc(MFt0, "NC_GLOBAL", "start_time")
  }else if(is(WEL, "RNetCDF")){
    if(!"Wells" %chin% var.get.nc(WEL, "outvars"))
      stop("wellC: NetCDF data set given for WEL does not contain Wells output")

    if(missing(welref)) welref <- as.data.table(unique({
      which.nc(WEL, "Wells", `!=`, 0, arr.ind = TRUE)[, 1:3]
    }, MARGIN = 1L))
    setnames(welref, c("C", "R", "L"))
  }else stop("wellC: invalid WEL; should be WEL.MFpackage or NetCDF")

  welref <- extract.welref(welref, dr, "wellC")

  # MODFLOW timings
  # - if WEL is given as a WEL.MFpackage, then stress period timings are
  #    sought
  mft <- switch(class(WEL)[1L],
                WEL.MFpackage = mfsptime(DIS, TRUE, MFt0),
                NetCDF = mftstime(WEL, TRUE))

  J <- wellJ(dr, welref)

  Q <- rep(list(double(length(mft) - 1L)), nrow(welref))
  names(Q) <- welref[, paste0("C", C, "R", R, "L", L)]
  switch(class(WEL)[1L],
         WEL.MFpackage = WEL$data[welref, {
           Q[[paste0("C", C, "R", R, "L", L)]][sp] <<- -sum(Q)
         }, on = c("C", "R", "L"), by = c("sp", "C", "R", "L")],
         NetCDF = welref[, {
           Q[[paste0("C", C, "R", R, "L", L)]] <<-
             -var.get.nc(WEL, "Wells", c(C, R, L, NA), c(1L, 1L, 1L, NA))
         }, by = c("C", "R", "L")])

  # relate the DRW time divisions to the MODFLOW time divisions
  drt.mft <- cellref.loc(dr@time, mft)

  # C = J/Q
  # the correct value of Q must be chosen for the time
  conc <- Map(`/`, J, lapply(Q, `[`, drt.mft))

  # turn Inf to NA
  for(i in 1:length(conc)) conc[[i]][!is.finite(conc[[i]])] <- NA_real_
  names(conc) <- names(J)

  conc
}

#' @import data.table
#'
extract.welref <- function(welref, dr, callingfun){
  # extract well reference
  # - should be data.table with C, R, L
  welref <- switch(class(welref)[1L],
                   matrix = as.data.table(welref),
                   data.frame = as.data.table(welref),
                   welref)
  switch(class(welref)[1L],
         data.table = unique(welref[, list(C, R, L)],
                             by = c("C", "R", "L")),
         numeric = data.table(C = as.integer(welref[1L]),
                              R = as.integer(welref[2L]),
                              L = as.integer(welref[3L])),
         integer = data.table(C = welref[1L],
                              R = welref[2L],
                              L = welref[3L]),
         WEL.MFpackage = unique(welref$data[, list(C, R, L)],
                                by = c("C", "R", "L")),
         RIV.MFpackage = unique(welref$data[, list(C, R, L)],
                                by = c("C", "R", "L")),
         GHB.MFpackage = unique(welref$data[, list(C, R, L)],
                                by = c("C", "R", "L")),
         DRN.MFpackage = unique(welref$data[, list(C, R, L)],
                                by = c("C", "R", "L")),
         STR.MFpackage = unique(welref$data[, list(C, R, L)],
                                by = c("C", "R", "L")),
         `NULL` = unique(dr@fluxout[J_out != 0, list(C, R, L)],
                         by = c("C", "R", "L")),
         stop(sprintf("%s: invalid welref", callingfun)))
}

#' @rdname DRWbalance
#' @name DRWbalance
#' @title Mass or Flux balance for a DRW model result
#'
#' @param dr
#' DRWmodel S4 class object;
#' results from \code{\link{DRW}}
#' @param type
#' \code{"mass"} or \code{"flux"}
#' @param plot
#' logical [1];
#' whether to plot the results using plot.DRWbalance before returning
#' @param balance;
#' a DRWbalance object, the result of \code{DRWbalance}
#' @param t
#' \code{NULL}, numeric [nts - 1] or a \code{\link{DRWmodel}};
#' time values to use for x axis; if \code{NULL} is given, then time step
#'  numbers are plotted
#' @param detail
#' logical [1];
#' whether to subdivide inputs, outputs and active mass
#' @param col
#' character [4];
#' colours to use for inputs, outputs, active mass and imbalance
#'  respectively
#' @param ...
#' \code{DRWbalance}: additional parameters to pass to
#'  \code{plot.DRWbalance};
#' \code{plot.DRWbalance}: axis range and label and title parameters
#'  (although there are internal defaults if these aren't given)
#'
#' @return
#' \code{DRWbalance} returns a list of class "DRWbalance", with elements:\cr
#' \code{$in_} numeric [nts - 1, 1]: input mass or flux from sources by
#'  time step\cr
#' \code{$out} numeric [nts - 1, 8]: output mass or flux to sinks,
#'  degradation or out of the model active region by time step; these
#'  values are negative\cr
#' \code{$Dactive} numeric [nts - 1, 2]: change in or rate of change of
#'  plume and sorbed mass by time step\cr
#' \code{$active} (only if \code{type == "mass"}) numeric [nts, 2]: mass in
#'  plume and sorbed by time step\cr
#' \code{$imbalance} numeric [nts - 1]: unaccounted for mass; should be 0
#'  apart from machine imprecision
#'
#' Each element has informative column names (try \code{colnames(bal$out)}
#'  for example).  \code{nts} means the number of time steps in the DRW
#'  model.  Because time step 1 of a DRW model represents the initial state,
#'  fluxes and changes of mass are between the time steps, which is why the
#'  \code{in_}, \code{out} and \code{Dactive} matrices have one less row
#'  than the number of time steps.
#'
#' The result is returned invisibly if \code{plot == TRUE}.
#'
#' \code{plot.DRWbalance} returns \code{NULL}.
#'
NULL

#' @rdname DRWbalance
#' @import data.table
#' @importFrom methods is
#' @export
#'
DRWbalance <- function(dr, type = c("mass", "flux")[1L], plot = FALSE,
                       t = dr, ...){
  if(!(isS4(dr) && is(dr, "DRWmodel")))
    stop("DRW::DRWbalance: dr must be a DRWmodel S4 class, the result of DRW")

  if(length(type) != 1L || !type %chin% c("mass", "flux"))
    stop("DRW::DRWbalance: type must be either \"mass\", or \"flux\", and length 1")

  nts <- length(dr@time)
  Dt <- dr@time[-1L] - dr@time[-nts]

  # inputs
  # - from sources
  in.srcM <- mapply(function(t1, t2, dt){
    sum(vapply(dr@release$J, do.call, double(100L),
               list(seq(t1 + dt/200, t2 - dt/200, length.out = 100L)))*
          dt/100, na.rm = TRUE)
  }, dr@time[-length(dr@time)], dr@time[-1L], Dt)
  if(type == "flux") in.srcJ <- in.srcM/Dt
  stopifnot(NROW(in.srcM) == nts - 1L)

  # outputs
  # - lost mass: dr@lost
  #  -- includes mass lost from any edge, into inactive cells or degraded
  out.lostM <- -dr@lost[-1L,]
  if(type == "flux") out.lostJ <- out.lostM/Dt
  stopifnot(NROW(out.lostM) == nts - 1L)
  #
  # - sinks
  out.snkJ <- double(nts - 1L)
  dr@fluxout[, out.snkJ[ts - 1L] <<- -sum(J_out), by = ts]
  if(type == "mass") out.snkM <- out.snkJ*Dt
  stopifnot(NROW(out.snkJ) == nts - 1L)

  # plume and sorbed mass
  plM <- double(nts)
  soM <- double(nts)
  if(nrow(dr@plume)) dr@plume[, plM[ts] <<- sum(m), by = ts]
  if(nrow(dr@sorbed)) dr@sorbed[, soM[ts] <<- sum(m), by = ts]
  #
  # - change in mass
  DplM <- diff(plM)
  DsoM <- diff(soM)
  if(type == "flux"){
    DplJ <- DplM/Dt
    DsoJ <- DsoM/Dt
  }
  stopifnot(NROW(DplM) == nts - 1L)
  stopifnot(NROW(DsoM) == nts - 1L)

  bal <- list(in_ = cbind(source = switch(type,
                                          mass = in.srcM,
                                          flux = in.srcJ)),
              out = cbind(sink = switch(type,
                                        mass = out.snkM,
                                        flux = out.snkJ),
                          switch(type,
                                 mass = out.lostM,
                                 flux = out.lostJ)),
              Dactive = cbind(plume = switch(type,
                                             mass = DplM,
                                             flux = DplJ),
                              sorbed = switch(type,
                                              mass = DsoM,
                                              flux = DsoJ)),
              active = if(type == "mass"){
                cbind(plume = plM, sorbed = soM)
              })

  bal$imbalance <-
    rowSums(bal$in_) + rowSums(bal$out) - rowSums(bal$Dactive)

  class(bal) <- c(type, "DRWbalance")

  if(plot) plot.DRWbalance(bal, t, ...)

  if(plot) invisible(bal) else bal
}

#' @rdname DRWbalance
#' @import graphics
#' @export
#'
plot.DRWbalance <- function(balance, t = NULL, detail = TRUE,
                            col = c(in_ = "red", out = "blue",
                                    Dactive = "darkgreen",
                                    imbalance = "black"), ...){
  # get title parameters
  dots <- list(...)

  xlab <- if(!"xlab" %in% names(dots)){
    switch(class(t)[1L],
           `NULL` = "time step",
           "time")
  }else dots$xlab
  ylab <- if(!"ylab" %in% names(dots)){
    class(balance)[1L]
  }else dots$ylab
  main <- if(!"main" %in% names(dots)) "" else dots$main
  sub <- if(!"main" %in% names(dots)) "" else dots$sub

  # sort out colours
  if(is.null(names(col))){
    names(col) <- c("in_", "out", "Dactive", "imbalance")
  }
  if(!all(c("in_", "out", "Dactive", "imbalance") %in% names(col))){
    stop("DRW::plot.DRWbalance: invalid names for col (colours)")
  }

  # sort out time values
  t <- switch(class(t)[1L],
              `NULL` = seq(2L, nrow(balance$in_) + 1L),
              integer = as.numeric(t),
              numeric = t,
              DRWmodel = t@time[-1L],
              stop("DRW::plot.DRWbalance: invalid t"))
  if(length(t) != nrow(balance$in_)) stop({
    "DRW::plot.DRWbalance: t is incorrect length; should be number of time steps - 1 (or nrow(balance$in_))"
  })

  # determing plotting range
  maxval <- with(balance, max(max(in_), max(Dactive), max(imbalance)))
  minval <- with(balance, min(min(out), min(Dactive), min(imbalance)))

  xlim <- if(!"xlim" %in% names(dots)) range(t) else dots$xlim
  ylim <- if(!"ylim" %in% names(dots)) c(minval, maxval) else dots$ylim

  # plot and legend layout
  layout(matrix(2:1, 2L, 1L), heights = c(3, 1))
  on.exit(layout(t(1L)))
  opar <- par("mar")

  # legend
  par(mar = c(.1, .1, .1, .1))
  plot.new()
  if(detail){
    legend("center",
           legend = c("in: source",
                      "out: sink",
                      "out: degrade",
                      "out: inactive",
                      "out: edge/ other",
                      "change in plume mass",
                      "change in sorbed mass",
                      "imbalance"),
         col = col[c("in_", "out", "out", "out", "out",
                     "Dactive", "Dactive", "imbalance")],
         lty = c(1L, 1L, 2L, 4L, 3L, 1L, 2L, 1L),
         title = class(balance)[1L], ncol = 3L)
  }else{
    legend("center",
           legend = c("in", "out", "change in active mass", "imbalance"),
           col = col[c("in_", "out", "Dactive", "imbalance")], lty = 1L,
           title = class(balance)[1L], ncol = 2L)
  }

  # plot
  par(mar = opar)
  plot.new()
  plot.window(xlim, ylim)
  box(); axis(1L); axis(2L)
  title(xlab = xlab, ylab = ylab, main = main, sub = sub)
  abline(h = 0, col = "grey")

  if(detail){
    lines(t, balance$in_[, "source"], col = col["in_"], lty = 1L)
    lines(t, balance$out[, "sink"], col = col["out"], lty = 1L)
    lines(t, balance$out[, "degraded"], col = col["out"], lty = 2L)
    lines(t, balance$out[, "inactive"], col = col["out"], lty = 4L)
    lines(t,rowSums(balance$out[, c("front", "left", "back",
                                    "right", "other")]),
          col = col["out"], lty = 3L)
    lines(t, balance$Dactive[, "plume"], col = col["Dactive"], lty = 1L)
    lines(t, balance$Dactive[, "sorbed"], col = col["Dactive"], lty = 2L)
    lines(t, balance$imbalance, col = col["imbalance"], lty = 1L)
  }else{
    lines(t, rowSums(balance$in_), col = col["in_"])
    lines(t, rowSums(balance$out), col = col["out"])
    lines(t, rowSums(balance$Dactive), col = col["Dactive"])
    lines(t, balance$imbalance, col = col["imbalance"])
  }

  NULL
}
