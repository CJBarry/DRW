library(DRW)
context("plot")

test_that("plots state correctly", {
  library(data.table)

  disnm <- system.file("drw_mf_demo.dis", package = "DRW", mustWork = TRUE)
  basnm <- system.file("drw_mf_demo.bas", package = "DRW", mustWork = TRUE)
  welnm <- system.file("drw_mf_demo.wel", package = "DRW", mustWork = TRUE)

  dis <- Rflow::read.DIS(disnm)
  bas <- Rflow::read.BAS(basnm, dis)
  wel <- Rflow::read.WEL(welnm, dis)

  gccs <- Rflow::gccs(dis)
  grcs <- Rflow::grcs(dis)

  state <- data.table(ts = 3L,
                      x = sample(seq(400, 600, 1), 50L),
                      y = sample(seq(400, 600, 1), 50L),
                      L = 1L, zo = .5, m = 1)
  rel <- data.table(ts = 3L, x = 500, y = 500, L = 1L, zo = .5, m = 10)

  expect_silent(plotDRWstate(state, rel, 3L, 1L, bas, wel, gccs, grcs))
  if(interactive()){
    expect_equal(readline("looks correct? (y/n)"), "y")
  }
  expect_silent(plotDRWstate(state, rel, 3L, 2L, bas, wel, gccs, grcs))
  if(interactive()){
    expect_equal(readline("looks correct? (y/n)"), "y")
  }
})
