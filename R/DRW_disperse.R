# DRW package - split and random walk dispersion

#' Title
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
  mpdt[, phi := if(sym){
    # using R's recycling feature (note that there will always be an even
    #  number of rows at this stage)
    rep(runif(.N/2, -pi, pi), each = 2L) + c(0, pi)
  }else runif(.N, -pi, pi)]
  #
  # - displacement
  #  -- dimensionless
  #   --- a random normal distribution with sd of 1
  #   --- if symmetrical, it is ensured that dispersion pairs have the same
  #        dimensionless displacement
  mpdt[, disp_ := if(sym){
    rep(rnorm(.N/2), each = 2L)
  }else rnorm(.N)]

  # apply displacement
  mpdt[, c("x", "y") := {
    # here, real-world means in real-world co-ordinates, i.e. not
    #  dimensionless

    # real-world dispersion coefficient, allowing for anisotropic
    #  dispersion (longitudinal, transverse)
    # - multiplied by the velocity for velocity-dependent dispersion (avv
    #    was set to 1 earlier if vdepD = FALSE)
    DC <- avv*
      (dt*(D[1L]*cos(phi - traj))^2L + (D[2L]*sin(phi - traj))^2L)^.5

    # real-world dispersion length
    disp <- DC*disp_

    list(x + disp*cos(phi),
         y + disp*sin(phi))
  }]

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
