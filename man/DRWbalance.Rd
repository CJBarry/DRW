% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse.R
\name{DRWbalance}
\alias{DRWbalance}
\alias{plot.DRWbalance}
\title{Mass or Flux balance for a DRW model result}
\usage{
DRWbalance(dr, type = c("mass", "flux")[1L], plot = FALSE, t = dr, ...)

\method{plot}{DRWbalance}(balance, t = NULL, detail = TRUE, col = c(in_ =
  "red", out = "blue", Dactive = "darkgreen", imbalance = "black"), ...)
}
\arguments{
\item{dr}{DRWmodel S4 class object;
results from \code{\link{DRW}}}

\item{type}{\code{"mass"} or \code{"flux"}}

\item{plot}{logical [1];
whether to plot the results using plot.DRWbalance before returning}

\item{t}{\code{NULL}, numeric [nts - 1] or a \code{\link{DRWmodel}};
time values to use for x axis; if \code{NULL} is given, then time step
 numbers are plotted}

\item{...}{\code{DRWbalance}: additional parameters to pass to
 \code{plot.DRWbalance};
\code{plot.DRWbalance}: axis range and label and title parameters
 (although there are internal defaults if these aren't given)}

\item{detail}{logical [1];
whether to subdivide inputs, outputs and active mass}

\item{col}{character [4];
colours to use for inputs, outputs, active mass and imbalance
 respectively}

\item{balance;}{a DRWbalance object, the result of \code{DRWbalance}}
}
\value{
\code{DRWbalance} returns a list of class "DRWbalance", with elements:\cr
\code{$in_} numeric [nts - 1, 1]: input mass or flux from sources by
 time step\cr
\code{$out} numeric [nts - 1, 8]: output mass or flux to sinks,
 degradation or out of the model active region by time step; these
 values are negative\cr
\code{$Dactive} numeric [nts - 1, 2]: change in or rate of change of
 plume and sorbed mass by time step\cr
\code{$active} (only if \code{type == "mass"}) numeric [nts, 2]: mass in
 plume and sorbed by time step\cr
\code{$imbalance} numeric [nts - 1]: unaccounted for mass; should be 0
 apart from machine imprecision

Each element has informative column names (try \code{colnames(bal$out)}
 for example).  \code{nts} means the number of time steps in the DRW
 model.  Because time step 1 of a DRW model represents the initial state,
 fluxes and changes of mass are between the time steps, which is why the
 \code{in_}, \code{out} and \code{Dactive} matrices have one less row
 than the number of time steps.

The result is returned invisibly if \code{plot == TRUE}.

\code{plot.DRWbalance} returns \code{NULL}.
}

