library(Rsubread)

dir.create('../../../output/mouse/subread-index',recursive = TRUE)
buildindex(basename = "../../../output/mouse/subread-index/GRCm39.genome",
           reference = "../../../data/annotation/mm39/GRCm39.genome.fa.gz",memory = 64000)
sessionInfo()
