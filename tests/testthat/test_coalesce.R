library(DRW)
context("coalesce")

test_that("required packages installed", {
  ip <- rownames(installed.packages())
  expect_true(all(c("coalesce", "Rflow", "stats", "data.table") %in% ip))
})

test_that("behaviour with particles in inactive cells", {
  od <- getwd()
  setwd(wd <- system.file(package = "DRW"))
  library(data.table)

  mfdata <- RNetCDF::open.nc(system.file("drw_mf_demo.nc",
                                         package = "DRW", mustWork = TRUE))
  wtop <- Rflow::get.wtop.nc(mfdata, "drw_mf_demo_wtop.nc")

  # one particle in an inactive cell
  state <- data.table(ts = 17L, x = 25, y = 25,
                      L = 1L, zo = .5, m = 5)

  expect_silent({
    statec <- coalesceDRW(state, c(20, 2), 1e-3, 100L, mfdata, wtop, 13L)
  })
  expect_equal(names(statec), c("state", "loss"))
  expect_equal(statec$loss, 5)
  expect_equal(nrow(statec$state), 0L)

  # one particle that is outside the model
  state <- data.table::data.table(ts = 17L, x = -5, y = -5,
                                  L = 1L, zo = .5, m = 5)

  expect_silent({
    statec <- coalesceDRW(state, c(20, 2), 1e-3, 100L, mfdata, wtop, 13L)
  })
  expect_equal(names(statec), c("state", "loss"))
  expect_equal(statec$loss, 5)
  expect_equal(nrow(statec$state), 0L)

  # one particle in an inactive cell, with other particles
  state <- data.table::data.table(ts = 17L,
                                  x = c(25, 750, 625),
                                  y = c(25, 500, 825),
                                  L = 1L, zo = c(.3, .5, .7), m = 5)

  expect_silent({
    statec <- coalesceDRW(state, c(20, 2), 1e-3, 100L, mfdata, wtop, 13L)
  })
  expect_equal(names(statec), c("state", "loss"))
  expect_equal(statec$loss, 5)
  expect_equal(sum(statec$state$m), 10)
  expect_equal(nrow(statec$state), 2L)

  setwd(od)
})

test_that("coalescing with z-offset", {
  od <- getwd()
  setwd(wd <- system.file(package = "DRW"))
  library(data.table)

  mfdata <- RNetCDF::open.nc(system.file("drw_mf_demo.nc",
                                         package = "DRW", mustWork = TRUE))
  wtop <- Rflow::get.wtop.nc(mfdata, "drw_mf_demo_wtop.nc")

  # shouldn't coalesce
  state <- data.table(ts = 17L,
                      x = 625, y = 825,
                      L = 1L, zo = c(.1, .9), m = 5)

  expect_silent({
    statec <- coalesceDRW(state, c(20, 20), 1e-3, 100L, mfdata, wtop, 13L)
  })
  expect_equal(names(statec), c("state", "loss"))
  expect_equal(nrow(statec$state), 2L)
  expect_equal(sum(statec$state$m), 10)

  # should coalesce
  state <- data.table(ts = 17L,
                      x = 625, y = 825,
                      L = 1L, zo = c(.4, .6), m = 5)

  expect_silent({
    statec <- coalesceDRW(state, c(20, 20), 1e-3, 100L, mfdata, wtop, 13L)
  })
  expect_equal(names(statec), c("state", "loss"))
  expect_equal(nrow(statec$state), 1L)
  expect_equal(sum(statec$state$m), 10)

  setwd(od)
})
