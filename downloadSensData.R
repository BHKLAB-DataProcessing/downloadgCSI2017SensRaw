library(PharmacoGx)
library(compareDrugScreens)

#####################################################
######### Prepare published molecular data ##########
#####################################################
    
options(stringsAsFactors = FALSE)


data("gcsi.genomics")
data(gcsi.genomics.mask)
data("ccle.gcsi.long")
data("gcsi.genomics.feature.info")
data("gcsi.line.info")

geneMap <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv")
data.types <- sapply(strsplit(colnames(gcsi.genomics), split = "\\."), function(x) return(x[[1]]))


rownames(gcsi.genomics) <- gcsi.line.info[rownames(gcsi.genomics),"CellLineName"]


gcsi.genomics[gcsi.genomics.mask] <- NA
gcsi.genomics.feature.info$GeneID <- gsub("GeneID:","",gcsi.genomics.feature.info$GeneID)	     
gcsi.genomics.feature.info <- gcsi.genomics.feature.info[-which(gcsi.genomics.feature.info$Symbol==""),]
		     
rnaseq <- gcsi.genomics[,data.types=="vsd"]
colnames(rnaseq) <- gsub("vsd.GeneID:", "", colnames(rnaseq))

##Check:
any(duplicated(colnames(rnaseq)))
rnaseq  <- t(rnaseq)

cnv <- gcsi.genomics[,data.types=="cn"]

colnames(cnv) <- gsub("cn.GeneID:", "", colnames(cnv))
any(duplicated(colnames(cnv)))
cnv <- t(cnv)
rownames(cnv) <- gcsi.genomics.feature.info$Symbol[match(rownames(cnv), gcsi.genomics.feature.info$GeneID)]		     
cnv <- cnv[-which(duplicated(rownames(cnv))),]
cnv <- cnv[na.omit(rownames(cnv)),]		     

loh <- gcsi.genomics[,data.types=="loh"]

colnames(loh) <- gsub("loh.GeneID:", "", colnames(loh))
any(duplicated(colnames(loh)))
loh <- t(loh)

mut <- gcsi.genomics[,data.types=="mut"]

colnames(mut) <- gsub("mut.GeneID:", "", colnames(mut))
any(duplicated(colnames(mut)))
mut <- t(mut)
rownames(mut) <- gcsi.genomics.feature.info$Symbol[match(rownames(mut), gcsi.genomics.feature.info$GeneID)]		     
	     

mutp <- gcsi.genomics[,data.types=="mutp"]
colnames(mutp) <- gsub("mutp.GeneID:", "", colnames(mutp))
any(duplicated(colnames(mutp)))
mutp <- t(mutp)


hot <- gcsi.genomics[,data.types=="hot"]
colnames(hot) <- gsub("hot.GeneID:", "", colnames(hot))
any(duplicated(colnames(hot)))
hot <- t(hot)

cellInfo <- gcsi.genomics[,data.types=="cln"]
colnames(cellInfo) <- gsub("cln.", "", colnames(cellInfo))

rownames(gcsi.genomics.feature.info) <- gsub("GeneID:", "", rownames(gcsi.genomics.feature.info))

molecInfo <- data.frame(cellid=colnames(rnaseq), row.names=colnames(rnaseq))
molecInfo <- cbind(molecInfo, batchid=NA)		     
		     
#rnaseq <- ExpressionSet(rnaseq)
#pData(rnaseq) <- molecInfo
#fData(rnaseq) <- gcsi.genomics.feature.info[rownames(rnaseq),]
#annotation(rnaseq) <- "rna"

cnv <- ExpressionSet(cnv)
pData(cnv) <- molecInfo
geneInfoM <- geneMap[match(rownames(cnv),geneMap[ , "gene_name"]), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(geneInfoM) <- rownames(cnv) 
geneInfoM <- geneInfoM[rownames(cnv),]
colnames(geneInfoM) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
fData(cnv) <- geneInfoM		     
annotation(cnv) <- "cnv"

#loh <- ExpressionSet(loh)
#pData(loh) <- molecInfo
#fData(loh) <- gcsi.genomics.feature.info[rownames(loh),]
		     
mut <- ExpressionSet(mut)
pData(mut) <- molecInfo
geneInfoM <- geneMap[match(rownames(mut),geneMap[ , "gene_name"]), c("gene_id", "EntrezGene.ID", "gene_name", "gene_biotype")]
rownames(geneInfoM) <- rownames(mut) 
geneInfoM <- geneInfoM[rownames(mut),]
colnames(geneInfoM) <- c("EnsemblGeneId", "EntrezGeneId", "Symbol", "GeneBioType")
fData(mut) <- geneInfoM			     
annotation(mut) <- "mutation"

#symbol <- gcsi.genomics.feature.info[rownames(mutp),]
#rownames(mutp) <- symbol$Symbol
#mutp <- ExpressionSet(mutp)
#pData(mutp) <- molecInfo
#fData(mutp) <- gcsi.genomics.feature.info[gcsi.genomics.feature.info$Symbol %in% rownames(mutp),]
#rownames(fData(mutp)) <- rownames(mutp)

#hot <- ExpressionSet(hot)
#pData(hot) <- molecInfo
#fData(hot) <- gcsi.genomics.feature.info[rownames(hot),]

#save(rnaseq, cnv, loh, mut, mutp, hot, cellInfo, file="/pfs/out/gCSI_molData.RData")
save(cnv,mut,cellInfo, file="/pfs/out/gCSI_molData.RData")

#####################################################
########### 2017 DATA (OLD) - SENSITIVITY ###########
#####################################################
    


data("ccle.gcsi.long")


library(data.table)

ccle.gcsi.long <- data.table(ccle.gcsi.long)
setkey(ccle.gcsi.long, CellLine, Drug)
gcsi.long <- ccle.gcsi.long[Group=="gCSI"]


gcsi.long[,Group:=NULL]

gcsi.long[,CellLine := gcsi.line.info[CellLine, "CellLineName"]]

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

## NB: there are more computed experiments than row, so we use the computed values to create the sens_info tbl. 
## differences seem to be from the sr786 cell line

## Prepare published IC50 values

data(gcsi.ic50)


gcsi.ic50 <- melt(gcsi.ic50, na.rm=TRUE)

colnames(gcsi.ic50) <- c("cellid","drugid", "IC50")

gcsi.ic50$cellid <- gcsi.line.info[gcsi.ic50$cellid, "CellLineName"]

experiments <- paste(gcsi.ic50[,"cellid"], gcsi.ic50[,"drugid"], sep="_")


#### Where did these come from?
experiments[which(!experiments %in% gcsi.long[,exp] )]

experiments[duplicated(experiments)] <- paste(experiments[duplicated(experiments)], "rep2", sep="_")

rownames(gcsi.ic50) <- experiments


sensitivity.info <- gcsi.ic50[,c("cellid", "drugid")]

# extra_rows <- rownames(sensitivity.info)[which(!rownames(sensitivity.info) %in% rownames(raw.sensitivity))]


data(gcsi.mv)
gcsi.mv <- melt(gcsi.mv, na.rm=TRUE)

colnames(gcsi.mv) <- c("cellid","drugid", "Mean Viability")
gcsi.mv$cellid <- gcsi.line.info[gcsi.mv$cellid, "CellLineName"]


experiments3 <- paste(gcsi.mv[,"cellid"], gcsi.mv[,"drugid"], sep="_")

#### Where did these come from?
experiments3[which(!experiments3 %in% gcsi.long[,exp] )]

experiments3[duplicated(experiments3)] <- paste(experiments3[duplicated(experiments3)], "rep2", sep="_")

rownames(gcsi.mv) <- experiments3

stopifnot(all(rownames(gcsi.mv) %in% rownames(sensitivity.info)))
stopifnot(all(rownames(raw.sensitivity) %in% rownames(sensitivity.info)))

gcsi.ic50 <- gcsi.ic50[rownames(sensitivity.info),]
gcsi.mv <- gcsi.mv[rownames(sensitivity.info),]

published.profiles <- data.frame(mean.viability_published = gcsi.mv[,"Mean Viability"], 
								 ic50_published = gcsi.ic50[,"IC50"])
rownames(published.profiles) <- rownames(sensitivity.info)

save(sensitivity.info, raw.sensitivity, gcsi.long, published.profiles, file="/pfs/out/raw.sensitivity.RData")

raw.sensitivity <- raw.sensitivity

raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gcsi2017_raw_sens_", i, ".rds"))

}
