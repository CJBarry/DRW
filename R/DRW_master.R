# DRW package - master function

DRW <- function(rootname, description, mfdir = ".",
                mfdata, wtop, dis, bas, wel, hds, cbb, cbf,
                source.term, STna.rm = FALSE,
                porosity,
                start.t, end.t, dt,
                D, Rf = 1, lambda = 0, decay.sorbed = FALSE,
                cd, mm, minnp = 100L, maxnp, Ndp = 2L,
                load.init = FALSE, init,
                Kregion = "auto", smd, dKcell, nKlpMFl = 1L,
                nc.to.mf = 1L, mfdata.split = FALSE,
                plot.state = TRUE,
                time.mismatch.tol = 1e-3){

  od <- getwd()
  setwd(mfdir)
  on.exit(setwd(od), add = TRUE)

  # don't allow a 0 or negative time.mismatch.tol
  time.mismatch.tol <- abs({
    if(abs(time.mismatch.tol) == 0) 1e-10 else time.mismatch.tol
  })

  # get MODFLOW inputs and outputs
  # make wtop
  # set up time steps
  # set up and preallocate
  # load initial state if requested
  # execute
  # - time step initial state
  # - source releases
  # - advect
  # - sinks
  # - disperse
  # - coalesce
  # - identify lost mass
  # - plot
  # z calculation
  # weighted kernel smooth

  # MODFLOW inputs and outputs ----
  #
  # --------------------------------------------------------------------- #
  # Connects to MF result NetCDFs and reads DIS, BAS and WEL (if given)
  #  input packages.
  # --------------------------------------------------------------------- #
  #
  # - open NetCDFs
  mfdatal <- lapply(mfdata, function(x){
    switch(class(x)[1L],
           character = open.nc(x),
           NetCDF = x,
           stop("DRW: invalid mfdata"))
  })
  on.exit(l_ply(mfdatal, close.nc), add = TRUE)
  ndset <- length(mfdatal)
  #
  # - check nc.to.mf (which relates each element of mfdata to a MODFLOW
  #    model number)
  if(!identical(length(nc.to.mf), length(ndset)))
    stop("DRW: nc.to.mf is length ", length(nc.to.mf),
         " but mfdata is length ", length(ndset),
         "; in most cases nc.to.mf should be 1:length(mfdata)")
  #
  # - gather time steps for each MODFLOW NetCDF and check that the end of
  #    one model's time period matches the start of the next model's
  #  -- time steps
  mftstl <- lapply(mfdatal, mftstime, absolute = TRUE)
  if(`&&`(!is.na(time.mismatch.tol),
          !isTRUE(all.equal(sapply(mftstl[-ndset], last),
                            sapply(mftstl[-1L], `[`, 1L),
                            tolerance = time.mismatch.tol)))) stop({
                              "DRW: end times of MODFLOW models do not match the start times of successive models"
                            })
  #
  # - MODFLOW model start times and final end time
  dsett <- c(lapply(mftstl, `[`, 1L), last(last(mftstl)), recursive = TRUE)
  #
  # - read MODFLOW package files
  disl <- lapply(dis, function(x) switch(class(x)[1L],
                                         character = read.DIS(x),
                                         DIS.MFpackage = x,
                                         stop("DRW: invalid dis")))
  basl <- Map(function(x, dis) switch(class(x)[1L],
                                      character = read.BAS(x, dis),
                                      BAS.MFpackage = x,
                                      stop("DRW: invalid bas")), bas, disl)
  well <- if(!missing(wel)) Map(function(x, dis){
    switch(class(x)[1L],
           character = read.WEL(x, dis),
           WEL.MFpackage = x,
           stop("DRW: invalid wel"))
  }, wel, disl)
  #
  # - gather stress period times for each MODFLOW model
  #  -- note that some elements may be repeated where there is not a simple
  #      one-to-one correspondence between NetCDF datasets and MODFLOW
  #      models, but this will not normally be a concern
  mfsptl <- Map(mfsptime, disl[nc.to.mf],
                lapply(mfdatal, att.get.nc, "NC_GLOBAL", "start_time"),
                MoreArgs = list(absolute = TRUE))
  #
  # - column and row co-ordinates (assumed to be the same for each MODFLOW
  #    model)
  gccs <- gccs(mfdata[[1L]], TRUE)
  grcs <- grcs(mfdata[[1L]], TRUE)


  # find the saturated groundwater top ----
  #
  # --------------------------------------------------------------------- #
  # see help(get.wtop.nc, Rflow)
  #
  # wtop may be a list of wtop NetCDFs or ready-made files; if wtop is
  #  not given, this code automatically assembles the correct wtop NetCDF
  # --------------------------------------------------------------------- #
  #
  wtopl <- if(missing(wtop)){
    Map(get.wtop.nc, mfdatal, paste0(rootname, "_wtop.nc"),
        MoreArgs = list(nts.dtit = if(mfdata.split) "sNTS" else "NTS"))
  }else{
    Map(function(m, w){
      switch(class(w),
             character = {
               get.wtop.nc(m, w, if(mfdata.split) "sNTS" else "NTS")
             },
             NetCDF = w,
             stop("DRW: invalid wtop"))
    }, mfdatal, wtop)
  }


  # set up DRW time steps ----
  #
  # --------------------------------------------------------------------- #
  # makes time steps with end time values such that:
  # 1. most time intervals are dt
  # 2. there are timesteps at the end of each mfdata dataset
  # 3. time values are ordered and not duplicated
  # 4. the time period for the DRW model does not extend beyond the MODFLOW
  #     models' time range
  # 5. the end time of the model is at most very slightly less than the end
  #     time of the last MODFLOW model, to avoid time step confusion
  # --------------------------------------------------------------------- #
  #
  tvals <- sort(unique(c(Map(seq, dsett[-ndset], dsett[-1L], dt),
                         start.t, dsett, end.t,
                         recursive = TRUE)))
  tvals <- tvals[tvals >= max(start.t, dsett[1L]) &
                   tvals <= min(end.t, last(dsett))]
  if(isTRUE(all.equal(last(tvals), last(dsett),
                      tolerance = time.mismatch.tol))){
    tvals[length(tvals)] <- tvals[length(tvals)] - time.mismatch.tol
  }
  ndrts <- length(tvals)
  if(tvals[1L] > start.t) warning({
    "DRW: simulation is starting after start.t because the MODFLOW models do not go far enough back"
  })
  if(last(tvals) < end.t - time.mismatch.tol*2) warning({
    "DRW: simulation is stopping before end.t because the MODFLOW models do not go far enough forward"
  })


  # load initial state from another model, if requested ----
  #
  # --------------------------------------------------------------------- #
  # Takes another DRWmodel and finds the time step just before tvals[1L],
  #  using this time step to give the plume (and sorbed) starting state for
  #  the current model.  Modifies tvals as necessary.
  # --------------------------------------------------------------------- #
  #
  if(load.init){
    # find which time step to use for initial state
    ts.init <- max({
      tmp <- which(init@time < tvals[1L])
      if(!length(tmp))
        stop("DRW: init starts too late to be used as the starting state")
      tmp
    })
    t.init <- init@time[ts.init]
    #
    # modify tvals or make error if init doesn't have a far enough advanced
    #  time step
    if(t.init < dsett[1L] - 4*time.mismatch.tol){
      stop("DRW: latest time in init@time is too far before the MODFLOW start time to be used as the initial state; 4*time.mismatch.tol before the MODFLOW start time is allowed")
    }else if(t.init >= dsett[1L] - 4*time.mismatch.tol &&
             t.init < dsett[1L] + 4*time.mismatch.tol){
      t.init <- dsett[1L]
    }
    tvals <- sort(unique(c(seq(t.init, tvals[1L], dt), tvals, dsett)))
    tvals <- tvals[tvals >= max(t.init, dsett[1L]) &
                     tvals <= min(end.t, last(dsett))]
    ndrts <- length(tvals)
  }


  # set up and preallocate ----
  #
  # --------------------------------------------------------------------- #
  # set up source term
  # preallocate lists and vectors to hold results
  # - mob
  # - immob
  # - fluxout
  # - degraded
  # - lost
  # --------------------------------------------------------------------- #
  #
  # - source term
  rel <- switch(class(source.term)[1L],
                data.table = source.term,
                DNAPLSourceTerm = ST.DNAPL(source.term),
                list = rbindlist(lapply(source.term, function(x){
                  switch(class(x)[1L],
                         data.table = x,
                         DNAPLSourceTerm = ST.DNAPL(x),
                         stop("DRW: one element of source.term is not valid"))
                })),
                stop("DRW: source.term is not valid"))
  #
  mob <- vector("list", ndrts)
  immob <- vector("list", ndrts)
  fluxout <- vector("list", ndrts)
  degraded <- numeric(ndrts)
  lost <- numeric(ndrts)
  #
  # - initial state
  mob[[1L]] <- if(load.init) init@plume[ts == ts.init]
  immob[[1L]] <- if(load.init && is.data.table(init@sorbed))
    init@sorbed[ts == ts.init]


  # execute ----
  #
  # --------------------------------------------------------------------- #
  # for each time step:
  # 1. bring forward state at end of last time step
  # 2. add source releases
  # 3. sorb and desorb
  # 4. advect (MODPATH)
  # 5. calculate sink and reaction fluxes (MassTrack)
  # 6. disperse
  # 7. register lost mass
  # 8. coalesce
  # 9. plot
  # --------------------------------------------------------------------- #
  #
  for(drts in 2:ndrts){

  }

}

#' DRW model results formal class
#'
#' @slot time numeric;
#' time values of the time step divides, including start and end of model
#' @slot plume data.table with columns:\cr
#' \code{$ts} (int, key): DRW time step\cr
#' \code{$x, $y, $z} (num): 3D position\cr
#' \code{$L} (int): MODFLOW layer\cr
#' \code{$zo} (num): z-offset within layer (0: cell base, to 1: top of
#'  groundwater in cell (see \code{\link[Rflow]{get.wtop.nc}}))\cr
#' \code{$m} (num): particle mass;
#' mobile particle swarm
#' @slot sorbed data.table, (columns as for \code{plume});
#' immobile particles
#' @slot release data.table with columns:\cr
#' \code{$x, $y, $L, $zo}: as with \code{plume}\cr
#' \code{$J} (list): list of functions giving flux in terms of time\cr
#' the source term release at each source location
#' @slot fluxout data.table with columns:\cr
#' \code{$ts} (int, key): DRW time step\cr
#' \code{$C, $R, $L} (int): MODFLOW cell reference\cr
#' \code{$J_out} (num): rate of mass loss in cell during this time step\cr
#' outflux of mass from particles by MODFLOW cell reference and DRW time
#'  step
#' @slot degraded numeric [ndrts];
#' mass lost to reactive degradation by DRW time step
#' @slot lost matrix [ndrts, 6];
#' mass lost from model by other means; each named column refers to mass
#'  lost from each edge of the MODFLOW model (top, left, right, bottom),
#'  failed source release due to release into an inactive or dry MODFLOW
#'  cell, or otherwise
#' @slot dispersion list with elements:\cr
#' \code{$D} named numeric [2 or 3]: longitudinal, transverse (and vertical)
#'  dispersivities (units of length) if \code{vdepD} or dispersion
#'  coefficients (units of length^2/time) otherwise\cr
#' \code{$vdepD} logical [1]: see above\cr
#' \code{$`3D`} logical [1]: \code{TRUE} if a 3D dispersion tensor is used
#' @slot reactions list with elements:\cr
#' \code{$Rf} numeric [1]: retardation factor\cr
#' \code{$lambda} numeric [1]: first-order reactive decay constant
#'  (equivalent to log(2)/half-life)\cr
#' \code{$decaysorbed} logical [1]: whether first-order reactive decay
#'  applies to immobile solute mass
#' @slot porosity numeric array [1], [NLAY] or [NCOL, NROW, NLAY];
#'  effective porosity of each MODFLOW cell, either uniform, by layer or by
#'  cell
#' @slot coalescing list with elements:\cr
#' \code{$cd} named numeric [2]: coalescing search radii horizontally and
#'  vertically\cr
#' \code{$mm} numeric [1]: minimum mass of particles after coalescing\cr
#' \code{$maxnp} integer [1]: maximumn number of particles after
#'  coalescing\cr
#' see (\code{\link[coalesce]{coalesce}}) for more details
#' @slot description character string;
#' description of the model run for future reference
#' @slot MFinfo character string [3];
#' descriptive attributes, including date, of the MODFLOW NetCDF, for future
#'  reference, so that it is clear what flow field (and what version) was
#'  used in the model
#' @slot run.timings POSIXct [];
#' summary of simulation timings, to give date of model run and for
#'  performance analysis
#'
#' @return
#' DRWmodel S4 object
#'
#' @export
#'
#' @examples
DRWmodel <- setClass("DRWmodel",
                     slots = c(time = "numeric",
                               plume = "data.table",
                               sorbed = "data.table",
                               release = "data.table",
                               fluxout = "data.table",
                               degraded = "numeric",
                               lost = "matrix",
                               dispersion = "list",
                               reactions = "list",
                               porosity = "array",
                               coalescing = "list",
                               description = "character",
                               MFinfo = "character",
                               run.timings = "POSIXct"))
