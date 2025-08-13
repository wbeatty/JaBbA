args<-commandArgs(TRUE)


library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(dryclean)


pon_object = pon$new(pon_path = "/gpfs/data/drazer-lab/Will_Jabba/resources/fixed.detergent.rds", wgs=TRUE)
dryclean_object <- dryclean$new(pon = pon_object)


fragcount_output <- args[1]
data <- readRDS(fragcount_output)
print(str(data))
print(summary(data))
clean <- dryclean_object$clean(cov=fragcount_output)
saveRDS(clean, args[2])
