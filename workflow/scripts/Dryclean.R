args<-commandArgs(TRUE)


library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(dryclean)




pon_object = pon$new(pon_path = "/gpfs/data/drazer-lab/Manoj_Jabba/fixed.detergent.rds", wgs=TRUE)
dryclean_object <- dryclean$new(pon = pon_object)


fragcount_output <- args[1]
clean <- dryclean_object$clean(cov=fragcount_output)
saveRDS(clean, args[2])

