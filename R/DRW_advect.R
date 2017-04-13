# DRW package - advection, using MODPATH 5

#' Particle Tracking Advection with MODPATH 5
#'
#' @param mpdt data.table; mobile particles
#' @param t1,t2 numeric [1]; time step start and end
#' @param MFx0,MFy0,MFt0 MODFLOW origin and start time
#' @param phi_e porosity
#' @param disnm file name for DIS package
#' @param dis DIS.MFpackage
#' @param bas BAS.MFpackage
#' @param hds file name for head save
#' @param cbb file name for cell-by-cell budget
#' @param cbf CBF file name
#' @param newdat write fresh DAT?
#' @param newcbf ask MODPATH to write fresh CBF?
#' @param transient does the MODFLOW model contain more than one time step?
#' @param maxnp number of particles MODPATH should allocate memory for
#'
#' @return
#' data.table with pathline results
#'
#' @import data.table
#'
advectMODPATH <- function(mpdt, t1, t2, MFx0, MFy0, MFt0, phi_e,
                          disnm, dis, bas, hds, cbb, cbf,
                          newds, newcbf, transient, maxnp){
  # write the DAT file if necessary
  # - if the model has reached a new MODFLOW dataset
  # - if the current DAT file doesn't allow for enough particles
  if(nrow(mpdt) > maxnp){
    write(dattxt(phi_e, nrow(mpdt)*2L, dis, bas), "DRW.dat")

    # update MPmaxnp back in the master function
    if(exists("MPmaxnp", parent.env(environment())))
      assign("MPmaxnp", nrow(mpdt)*2L, parent.env(environment()))
  }else if(newds) write(dattxt(phi_e, maxnp, dis, bas), "DRW.dat")

  # write MODPATH name file if necessary
  if(newds) write(namtxt(disnm, hds, cbb, "DRW.dat"), "DRW.nam")

  # no particles
  if(is.null(mpdt) || !nrow(mpdt)) return({
    data.table(ptlno = integer(0L), x = numeric(0L), y = numeric(0L),
               z_off = numeric(0L), z = numeric(0L), t = numeric(0L),
               C = integer(0L), R = integer(0L), L = integer(0L),
               timestep = integer(0L))
  })

  # write the response file
  write(rsptxt(t2 - MFt0, newcbf, cbf, transient), "DRW.rsp")

  # write starting locations file
  write(ptrtxt(mpdt[, list(x = x - MFx0, y = y - MFy0, L, zo)], t1 - MFt0),
        "DRW.ptr")

  # find MODPATH executable
  mpexe <- system.file("exec/Mpathr5_0.exe", package = "DRW",
                       mustWork = TRUE)

  # run MODPATH, noting execution time
  system(paste0(mpexe, " DRW.rsp"), intern = TRUE, wait = TRUE,
         show.output.on.console = FALSE)
  tm <- Sys.time()

  # read the pathline file, checking that it has been created recently
  if(!file.exists("pathline") || file.mtime("pathline") < tm - 5){
    stop("MODPATH failed")
  }
  ptl <- fread("pathline", skip = 1L)
  setnames(ptl, c("ptlno", "x", "y", "z_off", "z",
                  "t", "C", "R", "L", "timestep"))
  setkey(ptl, ptlno)
  ptl[, c("x", "y", "t") := list(x + MFx0, y + MFy0, t + MFt0)]
  ptl
}

# response file
arp <- "@RESPONSE\n"
rsptxt <- function(tlim, newcbf, cbf, transient){
  # intro line not needed
  paste0(arp, c("DRW.nam", # name file giving model data
                if(transient) "2", #
                if(transient) "1 0",
                "Y",
                paste(tlim, "1"), # terminate simulation at this time
                if(transient) ifelse(newcbf, "1", "2"), # new cbf?
                if(transient) cbf, # name of cbf to be written or read
                "2", # pathline output
                "N", # don't calculate location at specific time points
                "1",
                "DRW.ptr", # starting locations file
                rep("1", 2L),
                rep("N", 3L),
                "Y"), collapse = "\n")
}

#' MODPATH DAT file write
#'
#' @param por numeric; porosity
#' @param MXP integer; particles for which memory should be allocated
#' @param dis DIS.MFpackage
#' @param bas BAS.MFpackage
#'
#' @return
#' character string
#'
dattxt <- function(por, MXP, dis, bas){
  # initialise
  txt <- character(7L)

  txt[1] <- paste0(formatC(2^35L, 0L, 16L, "f"),
                   Rflow:::FFe(999, 16L, 6L, 3L),
                   Rflow:::FFe(1e30, 16L, 6L, 3L), "  ",
                   as.integer(MXP), "  1  1")

  txt[3] <- paste(dis$LAYCBD, collapse = " ")

  txt[4] <- paste(if(is.vector(ib <- bas$IBOUND)){
    vapply(ib, function(val) Rflow:::RIARRAY(CNSTNT = val, FMTIN_type = "i",
                                             FMTIN_w = 3L, flag.no = 10L),
           character(1))
  }else if(length(dim(ib)) == 2L){
    Rflow:::RIARRAY(arr = ib, FMTIN_type = "i", FMTIN_w = 3L, flag.no = 10L)
  }else{
    apply(ib, 3L, Rflow:::RIARRAY, FMTIN_type = "i",
          FMTIN_w = 3L, flag.no = 10L)
  }, collapse = "\n")

  # porosity is given either as a single value, a vector (value per layer)
  #  or an array (value per cell)
  if(is.vector(por)){
    lpor <- double(dis$extent["NLAY"])

    # fills layer values - works for uniform or value per layer
    lpor[] <- por
    txt[5] <- paste(vapply(lpor, function(p) Rflow:::RIARRAY(CNSTNT = p),
                           character(1L)), collapse = "\n")
  }else{
    txt[5] <- if(identical(length(dim(por)), 2L)){
      Rflow:::RIARRAY(arr = por, flag.no = 10L)
    }else{
      paste(apply(por, 3L, Rflow:::RIARRAY, flag.no = 10L), collapse = "\n")
    }
  }

  # TBEGIN - doesn't affect the running of MODPATH; this is the reference
  #  time for the start of the MODFLOW model
  txt[6] <- " 0.0"

  # starting and ending timesteps to process - it may be that using only the
  #  necessary timesteps will give a good saving in MODPATH run time; if you
  #  wish to develop this, then a new DAT file will be required for each DRW
  #  step
  txt[7] <- paste(c("", "1", "1", dis$extent["NPER"],
                    dis$sps[dis$extent["NPER"], "NSTP"]), collapse = "  ")

  return(paste(txt, collapse = "\n"))
}

#' MODPATH particle starting locations
#'
#' @param ptr.dat data.table (x, y, L, zo)
#' @param rt numeroc [1]; release time
#'
#' @return
#' character string
#'
#' @import data.table
#'
ptrtxt <- function(ptr.dat, rt){
  set(ptr.dat, NULL,
      c("C", "R", "It", "Jt", "Kt", "rt"),
      list(0L, 0L, 2L, 2L, 0L, rt))

  setcolorder(ptr.dat,
              c("C", "R", "L", "x", "y", "zo", "It", "Jt", "Kt", "rt"))

  ffmtptr <- mapply(formatC, ptr.dat,
                    width = c(4L, 4L, 3L, 17L, 17L, 17L, 2L, 2L, 2L, 13L),
                    digits = c(0L, 0L, 0L, 8L, 8L, 8L, 0L, 0L, 0L, 3L),
                    format = c("d", "d", "d", "e", "e",
                               "e", "d", "d", "d", "f"))

  if(nrow(ptr.dat) == 1L) paste(ffmtptr, collapse = "") else{
    apply(ffmtptr, 1L, paste, collapse = "")
  }
}

#' MODPATH name file text
#'
#' @param dis,hds,cbb,dat file names, not too long (use local names)
#'
#' @return
#' character string
#'
namtxt <- function(dis, hds, cbb, dat){
  paste(paste0("DIS    29    \'", dis, "\'"),
        paste0("main      10      \'", dat, "\'"),
        paste0("budget    17      \'", cbb, "\'"),
        paste0("head(binary)  18      \'", hds, "\'"), sep = "\n")
}
