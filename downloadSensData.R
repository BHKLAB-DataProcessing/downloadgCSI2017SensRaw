library(PharmacoGx)
# library(devtools)

# install.packages("BiocManager")
library(BiocManager)
message("trying install BiocManager")
message(.libPaths())
install(c("multtest", "NMF", "rmarkdown", "RColorBrewer"))

library(devtools)
message("trying install")
install_url("http://research-pub.gene.com/gCSI-cellline-data/compareDrugScreens_current.tar.gz")

# load("/pfs/gcsi2017raw/sensitivity.RData") #issue with Docker image. This is temp.

library(compareDrugScreens)

data("ccle.gcsi.long")


library(data.table)

ccle.gcsi.long <- data.table(ccle.gcsi.long)
setkey(ccle.gcsi.long, CellLine, Drug)
gcsi.long <- ccle.gcsi.long[Group=="gCSI"]


gcsi.long[,Group:=NULL]

gcsi.long[,exp := paste(gcsi.long[,CellLine], gcsi.long[,Drug], sep="_")]

areDuplicates <- duplicated(gcsi.long, by=c("exp", "Dose"))

gcsi.long[areDuplicates, exp:=paste(gcsi.long[areDuplicates, exp], "rep2", sep="_")]  

nconc <- gcsi.long[,.N,.(exp)]

max.con <- max(nconc[,N])

raw.sensitivity <- array(data = NA, dim=c(nrow(nconc), max.con, 2), dimnames=list(nconc[,exp],paste("Dose", 1:max.con), c("Dose", "Viability")))

setkey(gcsi.long, exp, Dose)

gcsi.long <- gcsi.long[order(exp, Dose)]

for (experiment in nconc[,exp]){
  
  raw.sensitivity[experiment,1:nconc[exp==experiment,N],] <- as.matrix(gcsi.long[exp==experiment,.(Dose, MedianViability)])
  
}

raw.sensitivity[,,"Viability"] <- raw.sensitivity[,,"Viability"]*100
gcsi.long[,`:=`(cellid = CellLine, drugid = Drug)]

sensitivity.info <- as.data.frame(gcsi.long[,.(unique(drugid),unique(cellid)), by=exp])
rownames(sensitivity.info) <- sensitivity.info$exp
colnames(sensitivity.info) <- c("expid", "drugid", "cellid")

save(sensitivity.info, raw.sensitivity, file="/pfs/out/raw.sensitivity.RData")

raw.sensitivity <- raw.sensitivity

raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gcsi2017_raw_sens_", i, ".rds"))

}
