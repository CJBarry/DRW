library(DRW)
context("master function")

test_that("basic functionality", {
  mfdir <- system.file(package = "DRW")
  od <- getwd()
  setwd(mfdir)
  mfdata <- RNetCDF::open.nc("drw_mf_demo.nc")

  # tests various forms for the source term (but not DNAPLSourceTerm)
  # - one of the point sources is into an inactive region
  source.term <- list({
    data.table::data.table(x = 625, y = 825,
                           L = 1L, zo = .5,
                           J = list(approxfun(c(0, 5000, 15000),
                                              c(0, 0, 1))))
  }, NULL, {
    # data frames don't accept list columns as easily as data tables
    df <- data.frame(x = 625, y = c(425, 975),
               L = 1L, zo = .5,
               J = 0)
    df$J <- list(approxfun(c(0, 5000, 15000),
                           c(0, 0, 1)))
    df
  })

  # tests file name and MFpackage forms for bas and wel input
  dis <- "drw_mf_demo.dis"
  bas <- switch(sample(1:2, 1L),
                "drw_mf_demo.bas",
                Rflow::read.BAS("drw_mf_demo.bas", dis))

  wel <- switch(sample(1:2, 1L),
                "drw_mf_demo.wel",
                Rflow::read.WEL("drw_mf_demo.wel", Rflow::read.DIS(dis)))

  expect_silent({
    testDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                   "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
                   "drw_mf_demo.dis", bas, wel,
                   "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                   newcbf = TRUE,
                   source.term = source.term,
                   porosity = matrix(.1, 30L, 20L),
                   start.t = 1000, end.t = 15000, dt = 1000,
                   D = c(10, 1), vdepD = TRUE,
                   cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                   Ndp = 4L)
  })

  # check that the saved object is the same as the returned object
  expect_equal(testDRW, readRDS("DRW_EXAMPLE.rds"))

  # time steps in this simple case
  expect_equal(testDRW@time, seq(1000, 15000, 1000), tolerance = .1)

  # test that mass loss to inactive cells is registered
  expect_gt((testDRW@lost)[length(testDRW@time), "inactive"], 900)
  expect_equal((testDRW@lost)[length(testDRW@time), "degraded"], 0,
               check.attributes = FALSE)

  # able to use old cbf
  expect_silent({
    testDRW2 <- DRW("DRW_EXAMPLE", "demo", mfdir,
                    "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
                    "drw_mf_demo.dis", bas, wel,
                    "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                    newcbf = FALSE,
                    source.term = source.term,
                    porosity = matrix(.1, 30L, 20L),
                    start.t = 1000, end.t = 15000, dt = 1000,
                    D = c(10, 1), vdepD = TRUE,
                    cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                    Ndp = 4L)
  })

  setwd(od)
})

test_that("reaction and dispersion processes", {
  mfdir <- system.file(package = "DRW")
  od <- getwd()
  setwd(mfdir)
  mfdata <- RNetCDF::open.nc("drw_mf_demo.nc")

  source.term <-
    data.table::data.table(x = 625, y = 825,
                           L = 1L, zo = .5,
                           J = list(approxfun(c(0, 5000, 15000),
                                              c(0, 0, 1))))

  dis <- "drw_mf_demo.dis"
  bas <- "drw_mf_demo.bas"
  wel <- "drw_mf_demo.wel"

  basicDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                  "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
                  "drw_mf_demo.dis", bas, wel,
                  "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                  newcbf = FALSE,
                  source.term = source.term,
                  porosity = matrix(.1, 30L, 20L),
                  start.t = 5000, end.t = 7500, dt = 100,
                  D = c(10, 1), vdepD = TRUE,
                  cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                  Ndp = 4L)

  # sorption
  expect_silent({
    testDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                   "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
                   "drw_mf_demo.dis", bas, wel,
                   "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                   newcbf = FALSE,
                   source.term = source.term,
                   porosity = matrix(.1, 30L, 20L),
                   start.t = 5000, end.t = 7500, dt = 100,
                   D = c(10, 1), vdepD = TRUE,
                   Rf = 10,
                   cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                   Ndp = 4L)
  })
  #
  # - check that the plume is less progressed than control
  expect_lt(testDRW@plume[ts == 10L, weighted.mean(x, m)],
            basicDRW@plume[ts == 10L, weighted.mean(x, m)])
  expect_gt(testDRW@plume[ts == 10L, weighted.mean(y, m)],
            basicDRW@plume[ts == 10L, weighted.mean(y, m)])

  # degradation
  expect_silent({
    testDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                   "drw_mf_demo.nc", "drw_mf_demo_wtop.nc",
                   "drw_mf_demo.dis", bas, wel,
                   "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                   newcbf = TRUE,
                   source.term = source.term,
                   porosity = matrix(.1, 30L, 20L),
                   start.t = 5000, end.t = 7500, dt = 100,
                   D = c(10, 1), vdepD = TRUE,
                   lambda = log(2)/1000,
                   cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                   Ndp = 4L)
  })
  #
  # - check that mass has been reduced

})
