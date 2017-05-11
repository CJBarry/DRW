# DRW package - split and random walk dispersion

#' Random Walk Dispersion
#'
#' @param mpdt
#' mobile particles data.table after the advection step
#' (pno, ts, x, y, L, zo, m, s, traj)
#' @param D
#' numeric [2 or 3];
#' simplified dispersion tensor (L, T, V)
#' @param vdepD
#' logical [1];
#' if \code{TRUE}, \code{D} represents the dispersivity (often denoted
#'  alpha, units of length), otherwise \code{D} represents the dispersion
#'  coefficient
#' @param dt
#' length of current time step
#' @param Ndp
#' integer [1];
#' number of dispersed particle pairs to spawn from each starting particle
#' @param sym
#' logical [1];
#' if \code{TRUE}, dispersed particle pairs are dispersed symmetrically,
#'  which ensures that centre of mass is absolutely conserved
#' @param ThreeDD
#' logical [1];
#' not yet coded
#'
#' @return
#' data.table (ts, x, y, L, zo, m)
#'
#' @import data.table
#' @importFrom pracma erfinv
#' @importFrom stats runif
#' @importFrom stats rnorm
#'
disperseRW <- function(mpdt, D, vdepD, dt, Ndp,
                       sym = TRUE, ThreeDD = FALSE){
  # 3D dispersion not yet programmed
  if(ThreeDD) warning("disperseRW: 3D dispersion not yet written")

  # calculate average velocity along each pathline
  # - if not vdepD, then D is taken as dispersion coefficient D, so avv is
  #    set to 1, so that it does not modify it
  mpdt[, avv := if(vdepD) s/dt else 1]

  # pathline trajectory should already be calculated from the pathline

  # expand table to accomodate new particles
  # - particles are replicated at this point
  mpdt <- mpdt[, .SD[rep(1L, 2L*Ndp)], by = pno]
  #
  # - divide mass between particles
  mpdt[, m := m/(2L*Ndp)]

  # assign random displacements and directions
  # - horizontal direction
  #  -- a random uniform distribution between -pi and pi radians
  #  -- if symmetrical, it is ensured that dispersion pairs are dispersed
  #      in exactly opposite directions
  mpdt[, xi2 := if(sym){
    # using R's recycling feature (note that there will always be an even
    #  number of rows at this stage)
    rep(runif(.N/2, -1, 1), each = 2L) + c(0, 1)
  }else runif(.N, -1, 1)]
  #
  # - displacement
  #  -- dimensionless
  #   --- a random normal distribution with sd of 1
  #   --- if symmetrical, it is ensured that dispersion pairs have the same
  #        dimensionless displacement
  mpdt[, xi1 := if(sym){
    rep(runif(.N/2, -1, 1), each = 2L)
  }else runif(.N, -1, 1)]

  mpdt[, rprime := erfinv(xi1)]
  mpdt[, phi := xi2*pi]
  mpdt[, c("Deltax", "Deltay") := {
    DL <- avv*D[1L]
    DT <- avv*D[2L]

    B <- 2^1.5*dt^.5
    dxl <- B*DL^.5*rprime*cos(phi)
    dxt <- B*DT^.5*rprime*sin(phi)

    # rotation applied
    list(dxl*cos(traj) - dxt*sin(traj),
         dxl*sin(traj) + dxt*cos(traj))
  }]
  mpdt[, c("x", "y") := list(x + Deltax, y + Deltay)]

  # return
  mpdt[, .(ts, x, y, L, zo, m)]
}

# this is not currently used
#' 2D rotation matrix
#'
#' @param phi numeric [1]; angle in radians
#' @param scale numeric [1]; expansion factor
#'
#' @return
#' 2 by 2 matrix
#'
Rot <- function(phi, scale = 1){
  scale*matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2L, 2L)
}
