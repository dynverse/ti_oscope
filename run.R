#!/usr/local/bin/Rscript

requireNamespace("dyncli", quietly = TRUE)
task <- dyncli::main()
# task <- dyncli::main(args = strsplit("--dataset ~/example.h5 --output ~/output.h5", " ")[[1]], definition_location = "ti_oscope/definition.yml")

library(dplyr, warn.conflicts = FALSE)
requireNamespace("dynutils", quietly = TRUE)
requireNamespace("dynwrap", quietly = TRUE)
# requireNamespace("Oscope", quietly = TRUE)
library(Oscope, warn.conflicts = FALSE)

#   ____________________________________________________________________________
#   Load data                                                               ####

counts <- t(as.matrix(task$counts))
parameters <- task$parameters

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####


# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = Sys.time())

# taken from vignette
# https://bioconductor.org/packages/release/bioc/vignettes/Oscope/inst/doc/Oscope_vignette.pdf
Sizes <- MedianNorm(counts+1, alternative = parameters$alternative_median)

DataNorm <- GetNormalizedMat(counts, Sizes)

MV <- CalcMV(
  Data = counts,
  Sizes = Sizes,
  MeanCutLow = parameters$mean_cut[[1]],
  MeanCutHigh = parameters$mean_cut[[2]],
  Plot = FALSE
)

DataSubset <- DataNorm[MV$GeneToUse,]

DataInput <- NormForSine(
  DataSubset,
  qt1 = parameters$qt[[1]],
  qt2 = parameters$qt[[2]]
)

SineRes <- OscopeSine(DataInput)

KMRes <- OscopeKM(
  SineRes,
  quan = parameters$quan
)

ToRM <- FlagCluster(SineRes, KMRes, DataInput)

KMResUse <- KMRes[-ToRM$FlagID]

if (length(KMResUse) == 0) {
  stop("No cyclic trajectory found")
}

ENIRes <- OscopeENI(
  KMRes = KMResUse,
  Data = DataInput,
  Ndg = parameters$ndg,
  NChun = parameters$nchun,
  N = parameters$niter,
  NCThre = parameters$ncthre
)

pseudotime <- percent_rank(order(ENIRes[[1]]))
names(pseudotime) <- colnames(counts)

# TIMING: done with method
checkpoints$method_aftermethod <- Sys.time()

#   ____________________________________________________________________________
#   Save output                                                             ####

output <-
  dynwrap::wrap_data(
    cell_ids = names(pseudotime)
  ) %>%
  dynwrap::add_cyclic_trajectory(
    pseudotime = pseudotime
  ) %>%
  dynwrap::add_timings(
    timings = checkpoints
  )

dyncli::write_output(output, task$output)
