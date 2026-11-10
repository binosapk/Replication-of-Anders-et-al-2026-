# Runs ONE numbered script of the v2 analysis:   Rscript run_step.R 04 solar     (output goes to stdout; redirect it)
# Step 00 builds state_base.RData. Every other step loads it and saves only the objects it creates
# to state_<step>_<sample>.RData, so several steps can run in parallel.
.a <- commandArgs(trailingOnly = TRUE); .step <- .a[1]; .samp <- if (length(.a) > 1) .a[2] else "all"
.V2 <- "C:/Users/bisap1788/OneDrive - Norwegian University of Life Sciences/Leeds/Solar-vs-Wind-Experiment/reanalysis_v2"
.ROOT <- dirname(.V2)
suppressPackageStartupMessages({ library(apollo); library(dplyr); library(tidyr) })
setwd(.V2)
if (.step != "00") {
  .e <- new.env(); load("state_base.RData", envir = .e)
  for (.nm in ls(.e, all.names = TRUE)) assign(.nm, get(.nm, .e), globalenv()); rm(.e, .nm)
  setwd(.ROOT)
}
SAMPLE <- .samp
.before <- ls(all.names = TRUE)
.file <- list.files(.V2, pattern = paste0("^", .step, "_.*\\.R$"), full.names = TRUE)
.t0 <- Sys.time()
.res <- try(source(.file, echo = FALSE, print.eval = TRUE, encoding = "UTF-8"))
.new <- setdiff(ls(all.names = TRUE), c(.before, ".res", ".new", ".out", ".t0"))
.out <- if (.step == "00") file.path(.V2, "state_base.RData") else file.path(.V2, sprintf("state_%s_%s.RData", .step, .samp))
if (.step == "00") save(list = setdiff(ls(all.names = TRUE), c(".a", ".step", ".samp", ".V2", ".ROOT", ".before", ".file", ".t0", ".res", ".new", ".out")), file = .out) else save(list = .new, file = .out)
cat(sprintf("\n[step %s (%s) finished in %.1f min]\n", .step, .samp, as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
if (inherits(.res, "try-error")) quit(status = 1)
