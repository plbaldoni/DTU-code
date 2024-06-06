OK <- requireNamespace(c('devtools','BiocParallel','doParallel'), quietly = TRUE)
if (!OK) {
  stop("devtools and BiocParallel packages required but is not installed (or can't be loaded)")
}
library(BiocParallel)
library(devtools)
library(doParallel)

load_all("../pkg")

workers <- 12

registerDoParallel(workers)
BPPARAM <- DoparParam()

summarizeSimulation(path = '../../output/simulation/data/',
                    dest = '../../output/simulation/summary/',
                    BPPARAM = BPPARAM)
