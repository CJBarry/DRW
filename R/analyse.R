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
  dr@fluxout[welref, {
    J[[paste0("C", C, "R", R, "L", L)]][ts] <<- sum(J_out)
  }, by = c("ts", "C", "R", "L"), on = c("C", "R", "L")]

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
