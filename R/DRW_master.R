# DRW package - master function

#' Dynamic Random Walk aqueous contaminant transport solution
#'
#' @param rootname
#' character string;
#' the rootname for the results file (\code{".rds"} is automatically
#'  affixed)
#' @param description
#' character string;
#' a description of this model run, for future reference
#' @param mfdir
#' character string;
#' path to the directory holding all the MODFLOW input and output files;
#'  the function navigates to this directory so that extended file paths
#'  for MODFLOW files needn't be specified
#' @param mfdata
#' NetCDF object, list of NetCDF objects or character string[Nmfds];
#' the MODFLOW data sets in NetCDF form (see \code{\link[Rflow]{GW.nc}}),
#'  either as NetCDF objects (see \code{\link[RNetCDF]{open.nc}}) or as
#'  character string file paths; this function may use a series of MODFLOW
#'  models with identical spatial co-ordinates and grids and with time
#'  periods that lead on from one another, thus a historic flow model may
#'  seamlessly be linked to a recent flow model, for example
#' @param wtop
#' NetCDF object, list of NetCDF objects or character string [Nmfds];
#' NetCDF data sets of the cell-by-cell top of water in each cell (the
#'  lower out of head and cell top); if character string file names are
#'  used, the files needn't exist beforehand as they can be created using
#'  \code{\link[Rflow]{get.wtop.nc}}
#' @param dis
#' character string [Nmfds];
#' file paths to DIS package files corresponding to the MODFLOW data sets
#'  in mfdata (note, this cannot be a list of DIS.MFpackage objects, as the
#'  file names are needed for MODPATH)
#' @param bas
#' BAS.MFpackage object (or list thereof) or character string [Nmfds];
#' the BAS packages corresponding to mfdata (see
#'  \code{\link[Rflow]{read.BAS}}); only the \code{$IBOUND} element is used
#' @param wel
#' Optional: WEL.MFpackage object (or list thereof) or character string
#'  [Nmfds];
#' the WEL packages corresponding to mfdata (see
#'  \code{\link[Rflow]{read.WEL}}); this information is only used for
#'  plotting, and plotting does not depend on it
#' @param hds,cbb
#' character string [Nmfds];
#' the head save and cell-by-cell budget files corresponding to mfdata;
#'  although this information is in mfdata, MODPATH needs to know where the
#'  original HDS files are
#' @param cbf
#' character string [Nmfds];
#' the composite budget file names created by MODPATH 5; if these files
#'  don't already exist or \code{newcbf = TRUE}, MODPATH 5 will write these
#'  files afresh
#' @param newcbf
#' logical [1];
#' Whether MODPATH 5 should rewrite the CBF file.  If you have modified the
#'  MODFLOW model since the last DRW run, then you should set this to
#'  \code{TRUE} (the default), but this can be time-consuming, so if you
#'  know it is not needed, set to \code{FALSE}.  The CBF will only ever be
#'  written for the first time step using a particular MODFLOW data set.
#' @param source.term
#' data.table, data.frame, \link[DNAPL]{DNAPLSourceTerm} object, or list of
#'  any combination of these;
#' information about the transient point releases in the system; if a data
#'  table or frame, columns should be: x, y (location in the same absolute
#'  co-ordinate system as \code{mfdata}), L (MODFLOW layer), zo (z-offset
#'  within layer) and J (list of functions of one variable, time, returning
#'  values representing source term flux, in units of mass per time OR
#'  numeric values representing constant release rates in units of mass per
#'  time)
#' @param STna.rm
#' logical [1];
#' if \code{TRUE}, any \code{NA} values resulting from the source term
#'  functions (within the model time frame) are ignored and treated as 0,
#'  otherwise an error will occur if \code{NA}s are found.
#' @param porosity
#' numeric [1], numeric [NLAY] or numeric array [NCOL, NROW, NLAY];
#' the porosity (fractional, 0 to 1), either as uniform value,
#'  layer-by-layer values or a 3D array of cell-by-cell values matching
#'  the MODFLOW grid dimensions
#' @param start.t,end.t
#' numeric [1];
#' start and end times of the model (the \code{\link[td]{td}} function may
#'  be useful)
#' @param dt
#' numeric [1]; time step size
#' @param D
#' numeric [2 or 3];
#' simplified dispersivity/ dispersion coefficient tensor, giving
#'  longitudinal, transverse and optionally vertical components
#' @param vdepD
#' logical [1];
#' \code{TRUE} if the dispersion coefficient is velocity-dependent, in
#'  which case \code{D} is taken to represent dispersivity, with units of
#'  length
#' @param Rf
#' numeric [1];
#' retardation factor (1 - (bulk_density times K_d)/(effective_porosity))
#' @param lambda
#' numeric [1]; first-order degradation constant (log(2)/half_life)
#' @param decay.sorbed
#' logical [1];
#' \code{TRUE} if first-order degradation affects sorbed contaminant (not
#'  yet implemented)
#' @param cd
#' numeric [2];
#' horizontal and vertical coalescing search radii (see
#'  \code{\link[coalesce]{coalesce}})
#' @param mm
#' numeric [1];
#' minimum mass of particles after coalescing (see
#'  \code{\link[coalesce]{coalesce}})
#' @param minnp
#' integer [1];
#' minimum number of particles on which to perform coalescing
#' @param maxnp
#' integer [1];
#' maximum number of particles after coalescing  (see
#'  \code{\link[coalesce]{coalesce}})
#' @param Ndp
#' integer [1];
#' number of pairs of particles to spawn to model the dispersion of
#'  contaminant mass
#' @param init
#' \code{NULL}, a \link{DRWmodel} object or a character string giving the
#'  RDS file name in which a DRWmodel object is stored;
#' another DRW model run to use as the initial state for the current model;
#'  the time step just before \code{start.t} is used as the initial state,
#'  modifying the time step times as necessary
#' @param nc.to.mf
#' integer [\code{length(mfdata)}];
#' the MODFLOW model that is related to by each NetCDF dataset in
#'  \code{mfdata}; normally this will be 1, 2, 3 ..., but sometimes a
#'  large MODFLOW model may be split into multiple NetCDF datasets (see
#'  \code{\link[Rflow]{GW.nc}}), in which case it will be something like
#'  1, 1, 1 ...; for the most common case in which there is only one
#'  MODFLOW model used, the default value of 1 is appropriate
#' @param plot.state
#' logical [1];
#' if \code{TRUE}, a summary plot of the model state will be shown after
#'  each time step (generally, not costly in terms of run time)
#' @param keep.MF.cellref
#' logical [2];
#' if \code{TRUE}, the final results for mobile and immobile particles will
#'  retain MODFLOW column (C), row (R), dataset (mfds) and time step (mfts)
#'  references
#' @param time.mismatch.tol
#' numeric [1];
#' various time comparisons are made during the model in order to check
#'  consistency, such as checking that the end time of one MODFLOW dataset
#'  matches the start time of the next; time.mismatch.tol gives the
#'  tolerance used in these comparisons in calls to
#'  \code{\link[base]{all.equal}}
#' @param ...
#' graphical parameters such as \code{xlim} (not axis labels, titles,
#'  \code{asp} or \code{zlim})
#' @param keepMPfiles
#' logical [1];
#' use for debugging when MODPATH not working: doesn't delete the MODPATH
#'  files
#'
#' @return
#' A \link{DRWmodel} object, invisibly.  The result is also saved to file:
#'  \code{paste0(rootname, ".rds")}.
#'
#' @import data.table
#' @import RNetCDF
#' @import Rflow
#' @import MassTrack
#' @import methods
#' @importFrom stats qunif
#' @export
#'
#' @examples
#' library(data.table)
#'
#' # single point source which ceases after t = 10000
#' mfdir <- system.file(package = "DRW")
#' demoDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
#'                "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
#'                "drw_mf_demo.dis", "drw_mf_demo.bas", "drw_mf_demo.wel",
#'                "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
#'                newcbf = TRUE,
#'                source.term = data.table(x = 625, y = 825,
#'                                         L = 1L, zo = .5,
#'                                         J = function(t){
#'                                           if(t < 1e4) 1 else 0
#'                                         }),
#'                porosity = .1,
#'                start.t = 9000, end.t = 11000, dt = 200,
#'                D = c(10, 1), vdepD = TRUE,
#'                cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
#'                Ndp = 4L)
#'
DRW <- function(rootname, description, mfdir = ".",
                mfdata, wtop, dis, bas, wel, hds, cbb, cbf, newcbf = TRUE,
                source.term, STna.rm = FALSE,
                porosity,
                start.t, end.t, dt,
                D, vdepD, Rf = 1, lambda = 0, decay.sorbed = FALSE,
                cd, mm, minnp = 100L, maxnp, Ndp = 2L,
                init = NULL,
                nc.to.mf = 1L,
                plot.state = TRUE,
                keep.MF.cellref = TRUE,
                time.mismatch.tol = 1e-3, ...,
                keepMPfiles = FALSE){

  run.start <- Sys.time()

  od <- getwd()
  setwd(mfdir)
  on.exit(setwd(od), add = TRUE)

  # set up an independent directory for MODPATH files
  # - this avoids confusion when multiple DRW models are running
  #    simultaneously
  # - the directory is based on rootname and is relative to mfdir, which
  #    has already been navigated to by this stage
  if(length(rootname) != 1L) stop({
    "DRW: `rootname` should be a length-1 character string (that is, one element, not one character)"
  })
  mpdir <- paste0(last(strsplit(rootname, "[\\/]")[[1L]]), "_DRW")
  if(!dir.exists(mpdir)) dir.create(paste0(mpdir))

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

  # MODFLOW inputs and outputs ----
  #
  # --------------------------------------------------------------------- #
  # Connects to MF result NetCDFs and reads DIS, BAS and WEL (if given)
  #  input packages.
  # --------------------------------------------------------------------- #
  #
  # - open NetCDFs
  if(is(mfdata, "NetCDF")) mfdata <- list(mfdata)
  #
  #  -- which elements of mfdata are given as character strings - these
  #      will need to be closed later
  mfdatawaschar <- sapply(mfdata, is.character)
  mfdatal <- lapply(mfdata, function(x){
    switch(class(x)[1L],
           character = open.nc(x),
           NetCDF = x,
           stop("DRW: invalid mfdata"))
  })
  on.exit(lapply(mfdatal[mfdatawaschar], close.nc), add = TRUE)
  ndset <- length(mfdatal)
  #
  # - check nc.to.mf (which relates each element of mfdata to a MODFLOW
  #    model number)
  if(!identical(length(nc.to.mf), ndset))
    stop("DRW: nc.to.mf is length ", length(nc.to.mf),
         " but mfdata is length ", ndset,
         "; in most cases nc.to.mf should be 1:length(mfdata)")
  #
  # - gather time steps for each MODFLOW NetCDF and check that the end of
  #    one model's time period matches the start of the next model's
  #  -- time steps
  mftstl <- lapply(mfdatal, mftstime, absolute = TRUE)
  if(`&&`(!is.na(time.mismatch.tol),
          !isTRUE(all.equal(sapply(mftstl[-ndset], last),
                            sapply(mftstl[-1L], `[`, 1L),
                            tolerance = time.mismatch.tol,
                            scale = 1)))) stop({
                              "DRW: end times of MODFLOW models do not match the start times of successive models"
                            })
  #
  # - MODFLOW model start times and final end time
  dsett <- c(lapply(mftstl, `[`, 1L), last(last(mftstl)), recursive = TRUE)
  #
  # - read MODFLOW package files
  #  -- the DIS files must be given as a character string file name,
  #      because the file name is required by MODPATH later
  disl <- lapply(dis, function(x) switch(class(x)[1L],
                                         character = read.DIS(x),
                                         stop({
                                           "DRW: invalid dis (must be character string file names)"
                                         })))
  basl <- if(is(bas, "BAS.MFpackage")) list(bas) else{
    Map(function(x, dis) switch(class(x)[1L],
                                character = read.BAS(x, dis),
                                BAS.MFpackage = x,
                                stop("DRW: invalid bas")), bas, disl)
  }
  well <- if(!missing(wel)) if(is(wel, "WEL.MFpackage")) list(wel) else{
    Map(function(x, dis){
      switch(class(x)[1L],
             character = read.WEL(x, dis),
             WEL.MFpackage = x,
             stop("DRW: invalid wel"))
    }, wel, disl)
  }
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
  gccs <- gccs(mfdatal[[1L]], TRUE)
  grcs <- grcs(mfdatal[[1L]], TRUE)
  #
  # - is each model a transient model (or, specifically, do they contain
  #    more than one time step)
  transientl <- sapply(disl,
                       function(x) sum(x$sps$NSTP) > 1L)
  #
  # - make sure that, if porosity is array, then it is 3D
  #  -- allows and corrects a 2D array if there is only one layer in the
  #      MODFLOW model
  if(is.array(porosity)){
    if(length(dim(porosity)) == 2L){
      if(dim.inq.nc(mfdatal[[1L]], "NLAY")$length != 1L){
        stop("DRW: porosity should be 3D if given as an array")
      }else porosity <- array(porosity, c(dim(porosity), 1L))
    }
    if(`||`(dim(porosity)[1L] != dim.inq.nc(mfdatal[[1L]], "NCOL")$length,
            dim(porosity)[2L] != dim.inq.nc(mfdatal[[1L]], "NROW")$length))
      stop("DRW: if porosity is an array, its dimensions should match the MODFLOW model [NCOL, NROW, NLAY]")
  }
  if(`&&`(is.vector(porosity),
          !(length(porosity) %in%
            c(1L, dim.inq.nc(mfdatal[[1L]], "NLAY")$length))))
    stop("DRW: if porosity is given as an R vector, then its length must be 1 (uniform value) or NLAY, (layer-by-layer values for each MODFLOW layer)")


  # find the saturated groundwater top ----
  #
  # --------------------------------------------------------------------- #
  # see help(get.wtop.nc, Rflow)
  #
  # wtop may be a list of wtop NetCDFs or ready-made files; if wtop is
  #  not given, this code automatically assembles the correct wtop NetCDF
  # --------------------------------------------------------------------- #
  #
  if(is(wtop, "NetCDF")) wtop <- list(wtop)
  # - determine if each MODFLOW NetCDF dataset is a subset and whether it
  #    is the first subset in a set
  mfdata.split <- vapply(mfdatal, function(nc){
    !is(try(att.inq.nc(nc, "NC_GLOBAL", "subset"), TRUE), "try-error")
  }, logical(1L))
  mfss1 <- vapply(mfdatal, function(nc){
    if(is(try(att.inq.nc(nc, "NC_GLOBAL", "subset"), TRUE),
          "try-error")){
      TRUE
    }else att.get.nc(nc, "NC_GLOBAL", "subset_start_ts") == 1L
  }, logical(1L))
  if(!mfss1[1L]) stop("DRW: the first subset of the MODFLOW set is not given, so the start time of the MODFLOW results cannot be determined")
  #
  # - determine the start times of the groundwater model sets
  #  -- this will be the same as dsett if none of the MODFLOW data sets are
  #      split
  #  -- if a data set is split, then the start time of the model set will
  #      be repeated for each subset
  mfsett <- numeric(ndset + 1L)
  for(i in 1:ndset) mfsett[i] <- {
    if(mfss1[i]) dsett[i] else mfsett[i - 1L]
  }
  mfsett[ndset + 1L] <- dsett[ndset + 1L]
  #
  wtopl <- if(missing(wtop)){
    wtopwaschar <- TRUE
    Map(get.wtop.nc, mfdatal, paste0(rootname, "_wtop.nc"),
        nts.dtit = ifelse(mfdata.split, "sNTS", "NTS"))
  }else{
    wtopwaschar <- sapply(wtop, is.character)
    Map(function(m, w, nd){
      switch(class(w),
             character = {
               get.wtop.nc(m, w, nd)
             },
             NetCDF = w,
             stop("DRW: invalid wtop"))
    }, mfdatal, wtop, ifelse(mfdata.split, "sNTS", "NTS"))
  }
  #
  on.exit(lapply(wtopl[wtopwaschar], close.nc), add = TRUE)


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
  tvals <- DRW_timesteps(start.t, end.t, dt, dsett, time.mismatch.tol)
  ndrts <- length(tvals)


  # load initial state from another model, if requested ----
  #
  # --------------------------------------------------------------------- #
  # Takes another DRWmodel and finds the time step just before tvals[1L],
  #  using this time step to give the plume (and sorbed) starting state for
  #  the current model.  Modifies tvals as necessary.
  # --------------------------------------------------------------------- #
  #
  load.init <- !is.null(init)
  #
  if(load.init){
    # read from file if necessary
    init <- switch(class(init)[1L],
                   DRWmodel = init,
                   character = readRDS(init),
                   stop("DRW: invalid init"))
    #
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
  # - lost
  # --------------------------------------------------------------------- #
  #
  # - source term
  rel <- switch(class(source.term)[1L],
                data.frame = as.data.table(source.term),
                data.table = source.term,
                DNAPLSourceTerm = ST.DNAPL(source.term, mfdatal, gccs, grcs),
                list = rbindlist(lapply(source.term, function(x){
                  switch(class(x)[1L],
                         data.frame = as.data.table(x),
                         data.table = x,
                         DNAPLSourceTerm = ST.DNAPL(x, mfdatal, gccs, grcs),
                         `NULL` = data.table(x = numeric(0L),
                                             y = numeric(0L),
                                             L = integer(0L),
                                             zo = numeric (0L),
                                             J = list()),
                         stop("DRW: one element of source.term is not valid"))
                })),
                `NULL` = data.table(x = numeric(0L), y = numeric(0L),
                                    L = integer(0L), zo = numeric (0L),
                                    J = list()),
                stop("DRW: source.term is not valid"))
  #
  #  -- allow constant value input for source.term$J
  if(is.numeric(rel$J)) rel$J <- lapply(rel$J,
                                        function(flux) function(t) flux)
  if(!all(c("x", "y", "L", "zo", "J") %in% names(rel)))
    stop("DRW: source.term does not have the correct columns")
  #
  mob <- vector("list", ndrts)
  immob <- vector("list", ndrts)
  fluxout <- vector("list", ndrts)
  lost <- matrix(0, ndrts, 7L,
                 dimnames = list(NULL,
                                 c("degraded", "inactive", "front", "left",
                                   "back", "right", "other")))
  #
  # - initial state
  if(load.init){
    mob[[1L]] <- init@plume[ts == ts.init]
    mob[[1L]][, ts := 1L]
    if(is.data.table(init@sorbed)){
      immob[[1L]] <- init@sorbed[ts == ts.init]
      immob[[1L]][, ts := 1L]
    }
  }
  #
  # - set the column orders to be used for the data tables
  #  -- during solution
  pcolorder <- c("ts", "x", "y", "L", "zo", "m")
  fcolorder <- c("ts", "C", "R", "L", "J_out")
  #
  # - only need a CBF file for a transient MODFLOW model (technically, one
  #    with more than one time step)
  newcbfl <- ifelse(transientl, newcbf, FALSE)
  #
  # - initial value
  mfds <- 0L
  #
  # - MODFLOW model origin
  MFx0l <- sapply(mfdatal, att.get.nc, "NC_GLOBAL", "origin-x")
  MFy0l <- sapply(mfdatal, att.get.nc, "NC_GLOBAL", "origin-y")
  if(uniqueN(MFx0l) != 1L || uniqueN(MFy0l) != 1L){
    stop("DRW: MODFLOW datasets have different origins")
  }
  MFx0 <- MFx0l[1L]; MFy0 <- MFy0l[1L]
  #
  # - an initial value to put in the MODPATH DAt file, allocating memory
  #    for pathlines
  MPmaxnp <- if(is.finite(maxnp)) maxnp*2L else 1e6L


  # execute ----
  #
  # --------------------------------------------------------------------- #
  # for each time step:
  #  1. bring forward state at end of last time step
  #  2. add source releases
  #  3. sorb and desorb
  #  4. advect (MODPATH)
  #  5. calculate sink and reaction fluxes (MassTrack)
  #  6. disperse
  #  7. register lost mass
  #  8. coalesce
  #  9. plot
  # 10. save
  # --------------------------------------------------------------------- #
  #
  tsattempt <- NULL
  for(drts in 2:ndrts) tsattempt <- try({
    if(is(tsattempt, "try-error")) break

    # time at start and end of stress period
    t1 <- tvals[drts - 1L]
    t2 <- tvals[drts]
    dift <- t2 - t1

    # establish MODFLOW model no., stress period and time step
    o.mfds <- mfds
    #
    # - start of time step
    mfds <- cellref.loc(t1, dsett)
    mfsp1 <- cellref.loc(t1, mfsptl[[mfds]])
    mfts1 <- cellref.loc(t1, mftstl[[mfds]])
    #
    # - end of time step
    #  -- NAs may be returned for values right at the end of a time
    #      division and these are corrected for
    mfsp2 <- cellref.loc(t2, mfsptl[[mfds]])
    if(is.na(mfsp2)) mfsp2 <- length(mfsptl[[mfds]]) - 1L
    #
    mfts2 <- cellref.loc(t2, mftstl[[mfds]])
    if(is.na(mfts2)) mfts2 <- length(mftstl[[mfds]]) - 1L

    # 1. bring state forward
    # - mobile
    statem <- copy(mob[[drts - 1L]])
    if(!is.null(statem)) statem[, ts := drts]
    #
    # - immobile
    statei <- copy(immob[[drts - 1L]])
    if(!is.null(statei)) statei[, ts := drts]

    # 2. source releases
    # - calculate released mass in this time step
    rel[, mtmp := vapply(J, function(f){
      sum(vapply(seq(t1 + dift/200, t2 - dift/200, length.out = 100L),
                 f, numeric(1L)), na.rm = STna.rm)*dift/100
    }, numeric(1L))]
    #
    # - error if NAs, but STna.rm will have wiped out any NAs if used
    if(any(is.na(rel$mtmp)))
      stop("DRW: time step ", drts, ": NAs in source release")
    #
    # - add released particles to current particle swarm
    statem <- rbind(if(is.data.table(statem)){
      statem[, list(ts, x, y, L, zo, m)]
    },
    rel[mtmp != 0, list(ts = drts,
                        x = x, y = y,
                        L = L, zo = zo,
                        m = mtmp)])
    #
    # - later development could allow release to the sorbed phase
    #
    # - remove unneeded columns for immobile phase
    #  -- very important the coalesce function needs to know to recalculate
    #      columns and rows after the sorbed material is added
    if(is.data.table(statei)){
      suppressWarnings(statei[, c("z", "C", "R", "mfds", "mfts") := NULL])
    }

    # 3. sorb and desorb
    if(Rf != 1){
      tmp <- sorb.desorb(copy(statem), copy(statei), Rf)
      statem <- tmp$mob
      statei <- tmp$immob
      rm(tmp)
    }

    # 4. advect
    # - write MODPATH DAT input?
    newds <- mfds != o.mfds
    #
    # - be safe: never read the wrong file
    if(newds){
      MPfiles <- c(paste0("DRW", c(".ptr", ".rsp", ".dat", ".nam")))
      file.remove(MPfiles[file.exists(MPfiles)])
    }
    #
    # - write MODPATH CBF input (composite budget file)?
    newcbf <- newds && newcbfl[nc.to.mf[mfds]]
    #
    # - MODFLOW start time
    MFt0 <- mfsett[mfds]
    #
    # - if this is a new MODFLOW data set, ensure that at least one
    #    particle is released in order to force the writing of the CBF file
    #    if required
    if(is.null(statem) || !nrow(statem)){
      statem <- data.table(ts = drts,
                           x = median(gccs), y = median(grcs),
                           L = 1L, zo = .5, m = 0)
    }
    #
    ptl <- advectMODPATH(copy(statem), t1, t2,
                         MFx0, MFy0, MFt0, porosity,
                         dis[nc.to.mf[mfds]],
                         disl[[nc.to.mf[mfds]]],
                         basl[[nc.to.mf[mfds]]],
                         hds[nc.to.mf[mfds]],
                         cbb[nc.to.mf[mfds]],
                         cbf[nc.to.mf[mfds]],
                         newds, newcbf,
                         transientl[[nc.to.mf[mfds]]],
                         MPmaxnp, mpdir)
    #
    # - correct time step number in case of split MODFLOW data set
    if(mfdata.split[mfds]){
      ptl[, timestep := timestep - as.integer({
        att.get.nc(mfdatal[[mfds]], "NC_GLOBAL", "subset_start_ts")
      }) + 1L]

      # When the end time of a MODPATH 5 simulation exactly corresponds
      #  with a MODFLOW time step divide, the particles are brought forward
      #  into the next time step and the z-offset adjusted for the head in
      #  the next time step.  With split MODFLOW NetCDFs, this causes an
      #  error with MassTrack which uses the NetCDF data set which won't
      #  then have the last time step in that NetCDF.  Similarly, when the
      #  start time is exactly a MODFLOW time step divide, then the previous
      #  time step is included with a trajectory of length and duration 0.
      #  This is not necessary and is simply deleted here (recognised as
      #  having time step 0 after the above correction).  MODPATH corrects
      #  the z-offset in this way when the head may change between time
      #  steps.
      if(any(ptl$timestep == 0L))
        ptl <- ptl[timestep != 0L]
      if(any(ptl$timestep >
             (snts <- dim.inq.nc(mfdatal[[mfds]], "sNTS")$length)))
        ptl <- ptl[timestep <= snts]
    }

    # 5. sinks and degradation
    if(nrow(ptl)){
      mt <- MassTrack(ptl, mfdatal[[mfds]], wtopl[[mfds]],
                      porosity, statem$m, TRUE, TRUE, FALSE,
                      TRUE, TRUE, FALSE, TRUE, lambda, TRUE, t2)
      #
      # register mass lost to sinks
      fluxout[[drts]] <- mt$traces[ml != 0,
                                   list(ts = drts,
                                        J_out = sum(ml)/dift),
                                   by = c("C", "R", "L")]
      #
      # register degraded mass
      lost[drts, "degraded"] <- sum(mt$traces$mrl)
      #
      # register mass that failed to release (inactive or dry cell)
      lost[drts, "inactive"] <- sum(mt$loss)
    }else mt <- list(traces = cbind(ptl, m = numeric(0L)))
    setnames(mt$traces, "z_off", "zo")
    #
    # - update statem, giving extra columns for pathline number, pathline
    #    length and average trajectory
    statem <- mt$traces[, c(list(pno = .GRP, ts = drts),
                            lapply(.SD, `[`, .N),
                            list(s = {
                              # sum of path segment lengths
                              sum((diff(x)^2L + diff(y)^2L)^.5)
                            }, traj = {
                              # trajectory from start to finish
                              atan2(y[.N] - y[1L], x[.N] - x[1L])
                            })),
                        by = ptlno, .SDcols = c("x", "y", "L", "zo", "m")]
    statem <- statem[m != 0]
    #
    # - decay sorbed phase if specified
    if(decay.sorbed) warning("DRW: sorbed phase degradation not yet programmed")

    # 6. disperse
    statem <- disperseRW(copy(statem), D, vdepD, dift, Ndp, TRUE, FALSE)

    # 7. register lost mass
    # - find column and row references
    statem[, C := cellref.loc(x, gccs)]
    statem[, R := cellref.loc(y, grcs, TRUE)]
    #
    # - register lost mass through each edge from particles with NA C or R
    #    references
    lost[drts, "front"] <- statem[is.na(R) & y <= grcs[1L], sum(m)]
    lost[drts, "left"] <- statem[is.na(C) & x <= gccs[1L], sum(m)]
    lost[drts, "back"] <- statem[is.na(R) & y >= last(grcs), sum(m)]
    lost[drts, "right"] <- statem[is.na(C) & x >= last(gccs), sum(m)]
    #
    # - remove lost particles
    statem <- statem[!is.na(C) & !is.na(R)]

    # 8. coalesce
    # - mobile
    if(!is.null(statem)){
      if(nrow(statem) > minnp){
        com <- coalesceDRW(statem, cd, mm, maxnp,
                           mfdatal[[mfds]], wtopl[[mfds]], mfts2)
        lost[drts, "inactive"] <- lost[drts, "inactive"] + com$loss
        statem <- com$state
        rm(com)
      }else{
        statem <- statem[, list(ts, x, y, L, zo, m)]
      }
    }
    #
    # - immobile
    if(!is.null(statei)){
      if(nrow(statei) > minnp){
        coi <- coalesceDRW(statei, cd, mm, maxnp,
                           mfdatal[[mfds]], wtopl[[mfds]], mfts2)
        lost[drts, "inactive"] <- lost[drts, "inactive"] + coi$loss
        statei <- coi$state
        rm(coi)
      }else{
        statei <- statei[, list(ts, x, y, L, zo, m)]
      }
    }

    # 9. plot if requested
    if(plot.state){
      # never stop running because this fails
      try(plotDRWstate(statem, rel, drts, mfsp2,
                       basl[[nc.to.mf[mfds]]], well[[nc.to.mf[mfds]]],
                       gccs, grcs, ...))
    }

    # 10. save
    if(!is.null(statem)) statem[, c("mfds", "mfts") := list(mfds, mfts2)]
    if(!is.null(statei)) statei[, c("mfds", "mfts") := list(mfds, mfts2)]
    #
    # - notation here avoids deleting elements if the mob or immob results
    #    are NULL
    mob[drts] <- list(statem)
    immob[drts] <- list(statei)
  })
  #
  # - unpack
  mob <- rbindlist(mob, use.names = TRUE, fill = TRUE)
  if(!length(mob)){
    mob <- data.table(ts = integer(0L), x = numeric(0L), y = numeric(0L),
                        L = integer(0L), zo = numeric(0L), m = numeric(0L),
                        mfds = integer(0L), mfts = integer(0L))
  }
  setkey(mob, ts)
  immob <- rbindlist(immob, use.names = TRUE, fill = TRUE)
  if(!length(immob)){
    immob <- data.table(ts = integer(0L), x = numeric(0L), y = numeric(0L),
                        L = integer(0L), zo = numeric(0L), m = numeric(0L),
                        mfds = integer(0L), mfts = integer(0L))
  }
  setkey(immob, ts)
  fluxout <- rbindlist(fluxout, use.names = TRUE, fill = TRUE)
  if(!length(fluxout)){
    fluxout <- data.table(ts = integer(0L),
                          C = integer(0L), R = integer(0L),
                          L = integer(0L), J_out = numeric(0L))
  }
  setkey(fluxout, ts)
  setcolorder(fluxout, fcolorder)
  #
  # - clean up (leave MODPATH summary file)
  if(!keepMPfiles){
    MPfiles <- paste0(mpdir, "/", c(paste0("DRW", c(".ptr", ".rsp", ".dat", ".nam")),
                                    "pathline", "endpoint"))
    file.remove(MPfiles[file.exists(MPfiles)])
  }


  # z calculation ----
  #
  # --------------------------------------------------------------------- #
  # determine z of each particle from layer and z-offset, given saturated
  #  thickness and elevation of each cell
  # also determines C and R references which are kept if keep.MF.cellref is
  #  set to TRUE (default)
  # --------------------------------------------------------------------- #
  #
  mob[, C := cellref.loc(x, gccs)]
  mob[, R := cellref.loc(y, grcs, TRUE)]
  immob[, C := cellref.loc(x, gccs)]
  immob[, R := cellref.loc(y, grcs, TRUE)]
  #
  # - get z
  if(nrow(mob)){
    mob[ts != 1L, z := {
      top.imtx <- cbind(C, R, L, mfts)
      bot.imtx <- cbind(C, R, L + 1L)

      top <- nc.imtx(wtopl[[mfds]], "wtop", top.imtx)
      bot <- nc.imtx(mfdatal[[mfds]], "elev", bot.imtx)

      HDRY <- att.get.nc(mfdatal[[mfds]], "Head", "HDRY")
      top[abs(top) < abs(HDRY)*1.00001 &
            abs(top) > abs(HDRY)*0.99999] <- NA_real_

      # qunif reads a quantile from a uniform distribution; it is fully
      #  vectorised
      qunif(zo, bot, top)
    }, by = mfds]
  }else mob[, z := numeric(0L)]
  if(nrow(immob)){
    immob[ts != 1L, z := {
      top.imtx <- cbind(C, R, L, mfts)
      bot.imtx <- cbind(C, R, L + 1L)

      top <- nc.imtx(wtopl[[mfds]], "wtop", top.imtx)
      bot <- nc.imtx(mfdatal[[mfds]], "elev", bot.imtx)

      HDRY <- att.get.nc(mfdatal[[mfds]], "Head", "HDRY")
      top[abs(top) < abs(HDRY)*1.00001 &
            abs(top) > abs(HDRY)*0.99999] <- NA_real_

      qunif(zo, bot, top)
    }, by = mfds]
  }else immob[, z := numeric(0L)]
  setcolorder(mob, c("ts", "x", "y", "z", "L", "zo", "m",
                     "C", "R", "mfds", "mfts"))
  setcolorder(immob, c("ts", "x", "y", "z", "L", "zo", "m",
                       "C", "R", "mfds", "mfts"))
  #
  # - remove MODFLOW cell reference information if it isn't wanted
  #  -- always keeps the layer reference as this is used for the Kernel
  #      smoothing post-processing operation
  if(!keep.MF.cellref){
    mob[, c("C", "R", "mfds", "mfts") := NULL]
    immob[, c("C", "R", "mfds", "mfts") := NULL]
  }


  run.end <- Sys.time()

  # save and return ----
  #
  # - get MODFLOW information
  MFinfo <- sapply(c("title", "author", "history"),
                   function(att){
                     sapply(mfdatal, att.get.nc,
                            "NC_GLOBAL", att)
                   })
  if(is.vector(MFinfo)) MFinfo <- t(MFinfo)
  #
  result <- DRWmodel(time = tvals,
                     plume = mob, sorbed = immob,
                     release = rel[, list(x, y, L, zo, J)],
                     fluxout = fluxout, lost = lost,
                     dispersion = mget(c("D", "vdepD", "Ndp")),
                     reactions = mget(c("Rf", "lambda", "decay.sorbed")),
                     porosity = as.array(porosity),
                     coalescing = mget(c("cd", "mm", "minnp", "maxnp")),
                     description = description,
                     MFinfo = MFinfo,
                     run.timings = c(run.start, run.end))
  #
  saveRDS(result, paste0(rootname, ".rds"))
  #
  invisible(result)
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
#' @slot lost matrix [ndrts, 7];
#' mass lost from model other than in sinks; each named column refers to
#'  mass lost from each edge of the MODFLOW model (top, left, right, bottom),
#'  failed source release due to release into an inactive or dry MODFLOW
#'  cell, degraded in first-order reactions, or other
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
#' @slot MFinfo character string matrix[number of MODFLOW models, 3];
#' descriptive attributes, including date, of the MODFLOW NetCDFs, for
#'  future reference, so that it is clear what flow field (and what version)
#'  was used in the model
#' @slot run.timings POSIXct [];
#' summary of simulation timings, to give date of model run and for
#'  performance analysis
#'
#' @return
#' DRWmodel S4 object
#'
#' @export
#'
DRWmodel <- setClass("DRWmodel",
                     slots = c(time = "numeric",
                               plume = "data.table",
                               sorbed = "data.table",
                               release = "data.table",
                               fluxout = "data.table",
                               lost = "matrix",
                               dispersion = "list",
                               reactions = "list",
                               porosity = "array",
                               coalescing = "list",
                               description = "character",
                               MFinfo = "matrix",
                               run.timings = "POSIXct"))
