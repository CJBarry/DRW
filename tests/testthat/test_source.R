library(DRW)
context("DNAPL source term")

if("DNAPL" %in% rownames(installed.packages())){
  library(DNAPL)

  test_that("DNAPL", {
    mfdir <- system.file(package = "DRW")
    od <- getwd()
    setwd(mfdir)
    mfdata <- RNetCDF::open.nc("drw_mf_demo.nc")

    cG <- cstG.DNmodel(.5, 20, 1, .5, .1, .1, .2, 1200, 1.2, 20, 3)

    fnm <- "DRW_EXAMPLE_dnst.rds"
    expect_silent(dnst <- {
      DNST(fnm, "test", cG,
           data.frame(year = c(1900, 1905), cons = c(0, 10000)),
           start.t = 0, end.t = 1500, dt = 10, x = 625, y = 825,
           mfdata = mfdata)
    })
    file.remove(fnm)

    dis <- "drw_mf_demo.dis"
    bas <- switch(sample(1:2, 1L),
                  "drw_mf_demo.bas",
                  Rflow::read.BAS("drw_mf_demo.bas", dis))

    wel <- switch(sample(1:2, 1L),
                  "drw_mf_demo.wel",
                  Rflow::read.WEL("drw_mf_demo.wel", Rflow::read.DIS(dis)))

    expect_silent({
      testDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                     mfdata, "drw_mf_demo_wtop.nc",
                     "drw_mf_demo.dis", bas, wel,
                     "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                     newcbf = TRUE,
                     source.term = dnst,
                     porosity = matrix(.1, 30L, 20L),
                     start.t = 1000, end.t = 15000, dt = 1000,
                     D = c(10, 1), vdepD = TRUE,
                     cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                     Ndp = 4L, STna.rm = TRUE)
    })

    # now with dnst@z0 as NA
    dnst@z0 <- NA_real_
    expect_silent({
      testDRW <- DRW("DRW_EXAMPLE", "demo", mfdir,
                     mfdata, "drw_mf_demo_wtop.nc",
                     "drw_mf_demo.dis", bas, wel,
                     "drw_mf_demo.hds", "drw_mf_demo.cbb", "DRWtest.cbf",
                     newcbf = TRUE,
                     source.term = dnst,
                     porosity = matrix(.1, 30L, 20L),
                     start.t = 1000, end.t = 15000, dt = 1000,
                     D = c(10, 1), vdepD = TRUE,
                     cd = c(20, 10), mm = 1e-7, minnp = 100L, maxnp = 2e4L,
                     Ndp = 4L, STna.rm = TRUE)
    })
  })
}
