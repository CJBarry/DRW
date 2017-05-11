library(DRW)
context("disperse")

test_that("required packages installed", {
  ip <- rownames(installed.packages())
  expect_true(all(c("stats", "data.table") %in% ip))
})

test_that("disperseRW", {
  for(trials in 1:5){
    # check one particle, visually
    state <- data.table::data.table(pno = 1L, ts = 8L,
                                    x = 0, y = 0,
                                    L = 1L, zo = .5, m = 10,
                                    s = stats::runif(1L, 2, 4),
                                    traj = stats::runif(1L, -pi, pi))

    D <- sample(3:24, 2L, TRUE)
    vdepD <- sample(c(TRUE, FALSE), 1L)
    Ndp <- 2000L
    sym <- sample(c(TRUE, FALSE), 1L)

    expect_silent(dstate <- disperseRW(state, D, vdepD, dt <- 2, Ndp, sym))
    expect_equal(names(dstate), c("ts", "x", "y", "L", "zo", "m"))

    # check number of particles
    expect_equal(nrow(dstate)/(2L*Ndp), nrow(state))

    # check mass conservation
    expect_equal(sum(dstate$m), sum(state$m), tolerance = 1e-10)

    # check centre of mass conservation if sym = TRUE
    if(sym) expect_equal(dstate[, {
      sapply(.SD, stats::weighted.mean, m)
    }, .SDcols = c("x", "y")], state[, {
      sapply(.SD, stats::weighted.mean, m)
    }, .SDcols = c("x", "y")])

    # visual checks
    if(interactive()){
      cat("\n")
      cat("D: "); print(D)
      cat("vdepD: "); print(vdepD)
      cat("Ndp: "); print(Ndp)
      cat("sym: "); print(sym)

      plot(dstate[, .(x, y)], col = "#FF000020", pch = 16L, asp = 1)
      lines(state[, list(x - c(s*cos(traj), 0),
                         y - c(s*sin(traj), 0))], lwd = 4)
      points(state[, .(x, y)], cex = 3, lwd = 4)
      legend("bottom", legend = c("start", "dispersed", "trajectory"),
             pt.cex = c(3, 1, NA), pt.lwd = c(1, 4, NA),
             lty = 1L, lwd = c(NA, NA, 4),
             col = c("black", "red", "black"), ncol = 3L)

      # plot one sd ellipse
      r0 <- ((if(vdepD) state$s/dt else 1)*D[1L]*2*dt)^.5
      r1 <- ((if(vdepD) state$s/dt else 1)*D[2L]*2*dt)^.5
      phi0 <- -state$traj
      lines(t(vapply(seq(-pi, pi, length.out = 101L),
                     function(phi){
                       c(matrix(c(cos(phi0), -sin(phi0),
                                  sin(phi0), cos(phi0)), 2L, 2L) %*%
                           c(r0*cos(phi), r1*sin(phi)))
                     }, numeric(2L))))

      expect_equal(readline("1sd ellipse looks correct? (y/n)"), "y")
    }
  }

  # check sd is roughly correct
  for(trial in 1:5){
    state <- data.table::data.table(pno = 1L, ts = 8L,
                                    x = 10, y = 5.6,
                                    L = 1L, zo = .5, m = 10,
                                    s = stats::runif(1L, 2, 4),
                                    traj = 0)

    D <- sample(3:24, 2L, TRUE)
    vdepD <- sample(c(TRUE, FALSE), 1L)
    Ndp <- 200000L
    sym <- sample(c(TRUE, FALSE), 1L)

    expect_silent(dstate <- disperseRW(state, D, vdepD, dt <- 2, Ndp, sym))

    expectsd <- ((if(vdepD) state$s/dt else 1)*D*2*dt)^.5

    expect_equal(c(sd(dstate$x), sd(dstate$y)), expectsd, tolerance = 1e-1)
  }
})
