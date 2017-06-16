# DRW package: time stepping

#' Set up time steps for DRW model
#'
#' @inheritParams DRW
#' @param dsett
#' numeric [];
#' MODFLOW dataset time divides, including start and end
#'
#' @import data.table
#' @return
#'
DRW_timesteps <- function(start.t, end.t, dt, dsett, time.mismatch.tol){

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
  ndset <- length(dsett) - 1L
  if(start.t > end.t) stop("DRW_steps: start time after end time")
  if(dt <= 0) stop("DRW_timesteps: dt <= 0")
  tvals <- sort(unique(c(Map(seq, dsett[-(ndset + 1L)], dsett[-1L], dt),
                         start.t, dsett, end.t,
                         recursive = TRUE)))
  tvals <- tvals[tvals >= max(start.t, dsett[1L]) &
                   tvals <= min(end.t, last(dsett))]
  if(isTRUE(all.equal(last(tvals), last(dsett),
                      tolerance = time.mismatch.tol, scale = 1))){
    tvals[length(tvals)] <- tvals[length(tvals)] - time.mismatch.tol
  }

  if(tvals[1L] > start.t) warning({
    "DRW_timesteps: simulation is starting after start.t because the MODFLOW models do not go far enough back"
  })
  if(last(tvals) < end.t - time.mismatch.tol*2) warning({
    "DRW_timesteps: simulation is stopping before end.t because the MODFLOW models do not go far enough forward"
  })

  if(any(is.na(tvals))) stop({
    "DRW_timestep: NA values in time steps"
  })
  if(is.unsorted(tvals, strictly = TRUE)) stop({
    "DRW_timestep: time steps are not in strictly ascending order"
  })

  tvals
}
