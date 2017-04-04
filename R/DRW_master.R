# DRW package - master function

DRW <- function(rootname, description, mfdir = ".",
                mfdata, wtop, dis, bas, wel, hds, cbb, cbf,
                source.term,
                porosity,
                start.t, end.t, dt,
                D, Rf = 1, lambda = 0, decay.sorbed = FALSE,
                cd, mm, maxnp, Ndp = 2L,
                load.init = FALSE, init,
                Kregion = "auto", smd, dKcell, nKlpMFl = 1L,
                nc.to.mf = 1L, mfdata.split = FALSE,
                plot.state = TRUE,
                time.mismatch.tol = 1e-3){

  od <- getwd()
  setwd(mfdir)
  on.exit(setwd(od), add = TRUE)

  # get MODFLOW inputs and outputs
  # make wtop
  # set up time steps
  # load initial state if requested
  # check source term
  # execute
  # - time step initial state
  # - source releases
  # - advect
  # - sinks
  # - disperse
  # - coalesce
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
    switch(class(x),
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
  disl <- lapply(dis, function(x) switch(class(x),
                                         character = read.DIS(x),
                                         DIS.MFpackage = x,
                                         stop("DRW: invalid dis")))
  basl <- Map(function(x, dis) switch(class(x),
                                      character = read.BAS(x, dis),
                                      BAS.MFpackage = x,
                                      stop("DRW: invalid bas")), bas, disl)
  well <- if(!missing(wel)) Map(function(x, dis){
    switch(class(x),
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
  tvals <- tvals[tvals >= start.t & tvals <= end.t]
  if(isTRUE(all.equal(last(tvals), last(dsett),
                      tolerance = time.mismatch.tol))){
    tvals[length(tvals)] <- tvals[length(tvals)] - dt/1e-3
  }
  ndrts <- length(tvals)
  if(tvals[1L] > start.t) warning({
    "DRW: simulation is starting after start.t because the MODFLOW models do not go far enough back"
  })
  if(last(tvals) < end.t) warning({
    "DRW: simulation is stopping before end.t because the MODFLOW models do not go far enough forward"
  })

}
