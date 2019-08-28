library(PharmacoGx)
library(devtools)

load("/pfs/gcsi2017raw/sensitivity.RData") #issue with Docker image. This is temp.


raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gcsi2017_raw_sens_", i, ".rds"))

}
