library(DRW)
context("advection")

test_that("required packages installed", {
  ip <- rownames(installed.packages())
  expect_true(all(c("Rflow", "MassTrack", "data.table") %in% ip))
})

test_that("MODPATH working", {
  od <- getwd()
  setwd(wd <- system.file(package = "DRW"))
  cat(wd, "\n")

  mpdir <- "tmp"
  if(!dir.exists(mpdir)) dir.create(mpdir)

  # the example MODFLOW model is part of the Rflow package
  hdsnm <- "drw_mf_demo.hds"
  cbbnm <- "drw_mf_demo.cbb"
  disnm <- "drw_mf_demo.dis"
  basnm <- "drw_mf_demo.bas"

  dis <- Rflow::read.DIS(disnm)
  bas <- Rflow::read.BAS(basnm, dis)

  # the mass is not used in this part of the solution, but should check
  #  that giving a data table with mass to advectMODPATH doesn't cause it
  #  to fail
  state <- data.table::data.table(x = 625, y = 825,
                                  L = 1L, zo = .5,
                                  m = 10)

  # the system command instructs not to print the MODPATH output, so this
  #  operation should be silent
  expect_silent({
    ptl <- DRW:::advectMODPATH(state, 1000, 2000, 0, 0, 0, .2, disnm,
                               dis, bas, hdsnm, cbbnm, "DRWtest.cbf",
                               TRUE, TRUE, TRUE, 100L, mpdir)
  })
  expect_equal(names(ptl), MassTrack:::PTL.headers[1:10])

  # now try again, using the existing DAT and CBF files
  expect_silent({
    ptl <- DRW:::advectMODPATH(state, 1000, 2000, 0, 0, 0, .2, disnm,
                               dis, bas, hdsnm, cbbnm, "DRWtest.cbf",
                               FALSE, FALSE, TRUE, 100L, mpdir)
  })
  expect_equal(names(ptl), MassTrack:::PTL.headers[1:10])

  # check expected behaviour with MFt0 and origin
  expect_silent({
    ptl <- DRW:::advectMODPATH(state, 1000, 2000, 0, 0, 1000, .2, disnm,
                               dis, bas, hdsnm, cbbnm, "DRWtest.cbf",
                               FALSE, FALSE, TRUE, 100L, mpdir)
  })
  expect_equal(names(ptl), MassTrack:::PTL.headers[1:10])
  #
  # - the particle was released at the (artificially imposed) MODFLOW start
  #    time in this case, so the first time value should be 1000 - 1000 = 0
  expect_equal(ptl$t[1L], 1000)
  #
  expect_silent({
    ptl <- DRW:::advectMODPATH(state, 1000, 2000, 10, 10, 0, .2, disnm,
                               dis, bas, hdsnm, cbbnm, "DRWtest.cbf",
                               FALSE, FALSE, TRUE, 100L, mpdir)
  })
  expect_equal(names(ptl), MassTrack:::PTL.headers[1:10])

  setwd(od)
})

test_that("release into inactive cell", {
  od <- getwd()
  setwd(wd <- system.file(package = "DRW"))

  mpdir <- "tmp"
  if(!dir.exists(mpdir)) dir.create(mpdir)

  # the example MODFLOW model is part of the Rflow package
  hdsnm <- "drw_mf_demo.hds"
  cbbnm <- "drw_mf_demo.cbb"
  disnm <- "drw_mf_demo.dis"
  basnm <- "drw_mf_demo.bas"

  dis <- Rflow::read.DIS(disnm)
  bas <- Rflow::read.BAS(basnm, dis)

  # (25, 25) is in a no-flow cell
  state <- data.table::data.table(x = c(25, 625), y = c(25, 825),
                                  L = 1L, zo = .5, m = 10)

  expect_silent({
    ptl <- DRW:::advectMODPATH(state, 1000, 2000, 0, 0, 0, .2, disnm,
                               dis, bas, hdsnm, cbbnm, "DRWtest.cbf",
                               FALSE, FALSE, TRUE, 100L, mpdir)
  })
  expect_equal(names(ptl), MassTrack:::PTL.headers[1:10])
  expect_equal(data.table::uniqueN(ptl$ptlno), 1L)

  setwd(od)
})
