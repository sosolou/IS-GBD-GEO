rm(list = ls()) 
geo.data.series.number = "GSE16561"
setwd("C:/Users/Administrator/Desktop/IS biomarkers/1.download")


Sys.setenv(LIB_XML = "$(MINGW_PREFIX)") 


library(Biobase)
library(GEOquery)
library(limma) 
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)

gset <- getGEO(geo.data.series.number, GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
samples <- colnames(exprs(gset))
geo.platform.number = gset@annotation    
probes_expr = exprs(gset);dim(probes_expr) 

phenoDat = pData(gset) 
phenoDat = phenoDat[c(40:63,1:39),]
write.csv(phenoDat, file=paste(geo.data.series.number,"_Type.csv",sep=""), row.names=F)

ex = exprs(gset)
library(impute)
imputed_gene_exp = impute.knn(ex,k=10,rowmax = 0.5,
                              colmax=0.8,maxp =3000, rng.seed=362436069)
ex = imputed_gene_exp$data
ex = as.data.frame(ex)
ex$ID = rownames(ex)
ids = as.data.frame(gset@featureData@data[["ID"]])
ids$Gene.symbol = gset@featureData@data[["Gene.symbol"]]  
colnames(ids) = c("ID","Gene.Symbol")

library(dplyr)
exprSet = inner_join(ids,ex,by = 'ID')         
exprSet = na.omit(exprSet)
exprSet = exprSet[-which(exprSet$Gene.Symbol==''),]

library(limma)
exprSet= avereps(exprSet[,-c(1,2)],             
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
exprSet = normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
exprSet = exprSet[,c(40:63,1:39)] 
library(tibble)
exprSet <- rownames_to_column(exprSet, var="geneNames")
exprSet <- exprSet[order(exprSet[,1]),]
write.table(exprSet, file=paste(geo.data.series.number,".csv",sep=""), row.names=F, col.names=T, quote=F, sep=",")


rm(list = ls()) 
geo.data.series.number = "GSE22255"
setwd("C:/Users/Administrator/Desktop/IS biomarkers/1.download")
Sys.setenv(LIB_XML = "$(MINGW_PREFIX)") 
library(Biobase)
library(GEOquery)
library(limma) 
library(gplots)
library(stringr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
gset <- getGEO(geo.data.series.number, GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
samples <- colnames(exprs(gset))
geo.platform.number = gset@annotation   
probes_expr = exprs(gset);dim(probes_expr)

phenoDat = pData(gset) 
write.csv(phenoDat, file=paste(geo.data.series.number,"_Type.csv",sep=""), row.names=F)
ex = exprs(gset)
library(impute)
imputed_gene_exp = impute.knn(ex,k=10,rowmax = 0.5,
                              colmax=0.8,maxp =3000, rng.seed=362436069)
ex = imputed_gene_exp$data
ex = as.data.frame(ex)
ex$ID = rownames(ex)
ids = as.data.frame(gset@featureData@data[["ID"]])
ids$Gene.symbol = gset@featureData@data[["Gene.symbol"]] 
colnames(ids) = c("ID","Gene.Symbol")
library(dplyr)
ids = ids[-grep('///',ids$Gene.Symbol),]        
exprSet = inner_join(ids,ex,by = 'ID')         
exprSet = na.omit(exprSet)
exprSet = exprSet[-which(exprSet$Gene.Symbol==''),]
library(limma)
exprSet= avereps(exprSet[,-c(1,2)],             
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
exprSet = normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
library(tibble)
exprSet <- rownames_to_column(exprSet, var="geneNames")
exprSet <- exprSet[order(exprSet[,1]),]
write.table(exprSet, file=paste(geo.data.series.number,".csv",sep=""), row.names=F, col.names=T, quote=F, sep=",")


rm(list = ls()) 
setwd("C:/Users/Administrator/Desktop/IS biomarkers/1.download")
library(tidyverse)
library(dplyr)
GSE1 <- read.table(file = "GSE22255.csv", header = TRUE, sep="," )
GSE2 <- read.table(file = "GSE16561.csv", header = TRUE, sep="," )
GSE <- inner_join(GSE1, GSE2, by = "geneNames")
GSE1_type <- read.table(file = "GSE22255_type.csv", header = TRUE, sep="," ) %>%
  dplyr::select(sample = geo_accession, gender = characteristics_ch1.1, group = "affected.status..disease.state..ch1") %>%
  mutate(data="GSE22255")

GSE2_type <- read.table(file = "GSE16561_type.csv", header = TRUE, sep="," ) %>%
  dplyr::select(sample = geo_accession, gender = characteristics_ch1, group = description) %>%
  mutate(data="GSE16561")


GSE_type <- rbind(GSE1_type,GSE2_type) %>%
  distinct()  %>%
  mutate(
    gender = case_when(
      gender == "gender: female" | gender == "gender: Female" ~ "female",
      gender == "gender: male" | gender == "gender: Male" ~ "male")
  ) %>%
  mutate(
    Group = case_when(
      group == "IS patient" | group == "Stroke" ~ "patient",
      group == "Control" | group == "control" ~ "normal")
  ) 

table(GSE_type$data,GSE_type$Group,GSE_type$gender)
library(sva)
library(limma)
rt <- as.matrix(GSE)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType <- c(rep(1,40),rep(2,63))
modType <- c(rep("Control",20),rep("IS",20),rep("Control",24),rep("IS",39))
mod <- model.matrix(~as.factor(modType))
normalize_GSE <- ComBat(data, batchType, mod, par.prior=TRUE)
normalize_GSE <- as.data.frame(normalize_GSE)
rm(mod,rt,batchType,data,dimnames,exp)
normalize_GSE_type <- GSE_type %>%
  arrange(Group, data) %>%
  mutate(id=paste0(sample,"_",Group))
normalize_GSE <- normalize_GSE[,c(normalize_GSE_type$sample)]
colnames(normalize_GSE) <- paste0(colnames(normalize_GSE),"_",normalize_GSE_type$Group) 
normalize_GSE_type_male <- normalize_GSE_type %>%
  filter(gender == "male")
normalize_GSE_type_female <- normalize_GSE_type %>%
  filter(gender== "female")
normalize_GSE_male <- normalize_GSE[,c((colnames(normalize_GSE) %in% normalize_GSE_type_male$id))]
normalize_GSE_female <- normalize_GSE[,(colnames(normalize_GSE) %in% normalize_GSE_type_female$id)]
normalize_GSE <- normalize_GSE %>% rownames_to_column("Gene") 
write.table(normalize_GSE, file="normalize_GSE.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type, file="normalize_GSE_type.txt", sep="\t", quote=F, col.names=T,row.names =F)
normalize_GSE_male <- normalize_GSE_male %>% rownames_to_column("Gene") 
write.table(normalize_GSE_male, file="normalize_GSE_male.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type_male, file="normalize_GSE_type_male.txt", sep="\t", quote=F, col.names=T,row.names =F)
normalize_GSE_female <- normalize_GSE_female %>% rownames_to_column("Gene") 
write.table(normalize_GSE_female, file="normalize_GSE_female.txt", sep="\t", quote=F, col.names=T,row.names =F)
write.table(normalize_GSE_type_female, file="normalize_GSE_type_female.txt", sep="\t", quote=F, col.names=T,row.names =F)


rm(list = ls())
library(tidyverse)
library(dplyr)
normalize_GSE_female <- read.table("normalize_GSE_female.txt", header=T, sep="\t", check.names=F) %>%
  column_to_rownames("Gene")
  
normalize_GSE_type_female <- read.table("normalize_GSE_type_female.txt", header=T, sep="\t", check.names=F)
table(normalize_GSE_type_female$data,normalize_GSE_type_female$Group)

conNum <- table(normalize_GSE_type_female$Group)[1]
treatNum <- table(normalize_GSE_type_female$Group)[2]


library(limma)
library(pheatmap)
normalize_GSE <- normalize_GSE_female

exp <- as.matrix(normalize_GSE) 
dimnames <- list(rownames(exp),colnames(exp))
rt <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt <- avereps(rt)

modType <- c(rep("normal",conNum),rep("patient",treatNum)) 
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("normal","patient")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(patient-normal,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff <- topTable(fit2,adjust='fdr',number=200000)
allDiffOut <- rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)


rm(design,cont.matrix,fit,fit2)


logFoldChange=0.58
adjustP=0.05

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSig <- diffSig[order(-diffSig$logFC),]
diffSigOut <- rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)


diffGeneExp=rt[rownames(diffSig),]
diffGeneExpOut=rbind(id=colnames(diffGeneExp),diffGeneExp)
write.table(diffGeneExpOut,file="diffGeneExp.txt",sep="\t",quote=F,col.names=F)


geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=normalize_GSE[hmGene,]
colnames(hmExp) = substr(colnames(hmExp),1,9)
Type=c(rep("Normal",conNum),rep("IS",treatNum))
names(Type)=substr(colnames(normalize_GSE),1,9)
Type=as.data.frame(Type)
pdf(file="heatmap_TOP50.pdf", width=10, height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


library(ggplot2)
Sig=ifelse((allDiff$adj.P.Val<adjustP) & (abs(allDiff$logFC)>logFoldChange), ifelse(allDiff$logFC>logFoldChange,"Up","Down"), "Stable")
rt=cbind(allDiff, Sig=Sig)

p <- ggplot(
  rt, aes(x = logFC, y = -log10(adj.P.Val), colour=Sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-logFoldChange,logFoldChange),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(adjustP),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15))

pdf(file="volcano.pdf", width=6, height=5.1)
print(p)
dev.off()

Incidence <- merge(Male,Female)
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)    
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

pvalueFilter=0.05       
qvalueFilter=0.05      

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


colnames(rt)[1]="Gene"
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)
showNum=8
if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="barplot_Top8.pdf", width=15, height=10)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()


pdf(file="bubble_Top8.pdf", width=20, height=15)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


ontology.col=c("#00AFBB", "#E7B800", "#90EE90")
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                  
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

rm(list = ls()) 
set.seed(123)
library(glmnet)                   
inputFile="C:/Users/Administrator/Desktop/IS biomarkers/2.diff/diffGeneExp.txt"       
inputFile1="C:/Users/Administrator/Desktop/IS biomarkers/2.diff/diff.txt"

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(data)

diff=read.table(inputFile1, header=T, sep="\t", check.names=F, row.names=1)
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit,cex.lab=1.5, cex.axis=1.5)
dev.off()
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LassoGene.txt", sep="\t", quote=F, row.names=F, col.names=F)
library(dplyr)
lassexpr <- data[c(lassoGene),] 
lassexprOut=rbind(id=colnames(lassexpr),lassexpr)
write.table(lassexprOut,file="LassoGeneExp.txt",sep="\t",quote=F,col.names=F)
library(limma)
library(pheatmap)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(lassexpr))
names(Type)=colnames(lassexpr)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=8, height=6)
pheatmap(lassexpr, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=12,
         fontsize_col=12)
dev.off()

Lasso.diff = diff[which(rownames(diff) %in% rownames(lassexpr)),]
data=t(lassexpr)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")
ddist=datadist(rt)
options(datadist="ddist")
library(rms)
library(rmda)
lrmModel=lrm(Type~ FOS+LCN2+F13A1+GPAT3+F5+RAB31+TNFSF13B+TNFSF10+FGL2+PDK4+SLC40A1+LAIR2+PLK2+TNFRSF4+FKBP11+MT1X+GPR18, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.5,0.99),
              lp=F, funlabel="Risk of Disease")
pdf("Nomo.pdf", width=10, height=8)
par(cex = 1.2)
plot(nomo)
dev.off()
par(cex = 1)
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()
rt$Type=ifelse(rt$Type=="normal", 0, 1)
dc=decision_curve(Type ~ FOS+LCN2+F13A1+GPAT3+F5+RAB31+TNFSF13B+TNFSF10+FGL2+PDK4+SLC40A1+LAIR2+PLK2+TNFRSF4+FKBP11+MT1X+GPR18, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()
lasso.prob <- predict(cvfit, newx=x , s=c(cvfit$lambda.min,cvfit$lambda.1se) )
y = ifelse(y=="normal", 0, 1)
re=cbind(y ,lasso.prob)


library(ROCR)
library(caret)
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")

tpr_min = performance(pred_min,"tpr")@y.values[[1]]
tpr_1se = performance(pred_1se,"tpr")@y.values[[1]]
dat = data.frame(tpr_min = perf_min@y.values[[1]],
                 fpr_min = perf_min@x.values[[1]]
)
library(ggplot2)
pdf(file="ROC.pdf", width=8, height=8)
ggplot() + 
  geom_line(data = dat,aes(x = fpr_min, y = tpr_min),color = "skyblue",size=1.5) + 
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "orange",size = 1.5,linetype=3 )+
  theme_bw()+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",round(auc_min,2)),color = "black",size=8)+
  scale_x_continuous(name  = "1-Specificity")+
  scale_y_continuous(name = "Sensitivity")+ 
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text=element_text(size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
dev.off() 



rm(list = ls()) 
library(randomForest)
library(dplyr)
set.seed(123)

inputFile="C:/Users/yuanfang/Desktop/ISbiomarkers/2.diff/female/diffGeneExp.txt"       
setwd("C:/Users/yuanfang/Desktop/ISbiomarkers/4.randomForest/female")    

inputFile1="C:/Users/yuanfang/Desktop/ISbiomarkers/2.diff/female/diff.txt"
diff=read.table(inputFile1, header=T, sep="\t", check.names=F, row.names=1)
rownames(diff)[which(rownames(diff)=="HLA-DQB1")] <- 'HLADQB1'
#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data) %>% as.data.frame() %>% rename(HLADQB1=`HLA-DQB1`)

data <- data[,which(colnames(data) %in% rownames(diff))]
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

data$Type=group

library(caret)
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]
train_group <- train$Type 
library(dplyr)
train1 <- train %>% 
  dplyr::select(-Type)
rf <- randomForest(as.factor(train_group)~., 
                   data=train1, 
                   ntree=500,
                   mtry=3,
                   importance=TRUE,
                   proximity=TRUE)
table(predict(rf),as.factor(train_group))

pdf(file="randomForest.pdf", width=6, height=6)
par(cex = 1.5)
plot(rf, main=" ", lwd=2)
dev.off()
par(cex = 1)

library(pROC)
pre = as.data.frame(predict(rf,test,type="prob"))
y=ifelse(test$Type=="normal", 0, 1)
roc1=roc(y, as.numeric(pre[,2]))
ciVec1=as.numeric(ci.auc(roc1, method="delong"))
pdf(file="ROC.pdf", width=6, height=6)
plot(smooth(roc1,method="density"), print.auc=F, legacy.axes=T, main="", col="#FF8884")
legend('bottom',
       c(paste0('   RF: ',sprintf("%.03f",roc1$auc),' (',paste0("95%CI: ",sprintf("%.03f",ciVec1[1]),"-",sprintf("%.03f",ciVec1[3])),')')
       ),
       col="#FF8884", lwd=2, bty = 'n')
dev.off()


auc_value <- auc(roc1)
print(auc_value)

optionTrees=which.min(rf$err.rate[,1])
optionTrees
write.table(optionTrees, file="optionTrees.txt", sep="\t", quote=F, col.names=F, row.names=F)

rf2=randomForest(as.factor(train_group)~., train1, ntree=optionTrees)

importance=importance(x=rf2)

pdf(file="geneImportance.pdf", width=8, height=8)
par(cex = 1.3)
varImpPlot(rf2, main="")
dev.off()
par(cex = 1)
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0.8])    
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)

library(limma)
library(pheatmap)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(sigExp))
names(Type)=colnames(sigExp)
Type=as.data.frame(Type)
sigExp=as.data.frame(sigExp)
pdf(file="heatmap.pdf", width=8, height=6)
pheatmap(sigExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =T,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=12,
         fontsize_col=12)
dev.off()

rf.diff = diff[which(rownames(diff) %in% rownames(sigExp)),]
paste(rownames(rf.diff), collapse="+")


data=t(sigExp)
rt=cbind(as.data.frame(data), Type=group)

library(rms)
ddist=datadist(rt)
options(datadist="ddist")

library(rms)
library(rmda)
lrmModel=lrm(Type~ PDK4+TNFSF10+ARHGEF40+TNFSF13B+CTSS, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.3,0.99),
              lp=F, funlabel="Risk of Disease")
pdf("Nomo.pdf", width=12, height=12)
par(cex = 1.3)
plot(nomo)
dev.off()
par(cex = 1)

cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

rt$Type=ifelse(rt$Type=="normal", 0, 1)
dc=decision_curve(Type ~ PDK4+TNFSF10+ARHGEF40+TNFSF13B+CTSS, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()


rm(list = ls()) 

library(venn)                   
setwd("C:/Users/yuanfang/Desktop/ISbiomarkers/5.intersect genes/female")    
geneList=list()

rt=read.table("C:/Users/yuanfang/Desktop/ISbiomarkers/3.lasso/female/LASSOGene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])       
uniqGene=unique(geneNames)        
geneList[["Lasso"]]=uniqGene

#读取SVM的结果文件
rt=read.table("C:/Users/yuanfang/Desktop/ISbiomarkers/4.randomforest/female/rfGenes.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])       
uniqGene=unique(geneNames)        
geneList[["Random Forest"]]=uniqGene


mycol=c("skyblue","pink")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()


intersectGene=Reduce(intersect,geneList)
paste(intersectGene,collapse=", ")
write.table(file="intersectGene.txt", intersectGene, sep="\t", quote=F, col.names=F, row.names=F)

diff=read.table("C:/Users/yuanfang/Desktop/ISbiomarkers/2.diff/female/diffGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)
intersectGeneExpr=diff[c(intersectGene),]
GeneExprOut=rbind(id=colnames(intersectGeneExpr),intersectGeneExpr)
write.table(file="intersectGeneExpr.txt", GeneExprOut, sep="\t", quote=F, col.names=F, row.names=T)


library(ggvenn) 
pdf(file="ggvenn.pdf", width=5, height=5)
ggvenn(geneList, show_elements = TRUE,fill_color = c("skyblue", "pink"),
       label_sep = "\n", stroke_size = 0.5,set_name_size = 4,
       text_size = 2.5)
dev.off()

rm(list = ls()) 

setwd("C:/Users/yuanfang/Desktop/ISbiomarkers/5.intersect genes/female")     

data=read.table("intersectGeneExpr.txt" , header=T, sep="\t", check.names=F, row.names=1) 

data=as.matrix(scale(t(data)))                                   
data.class <- rownames(data)
data.pca <- princomp(data,cor = T) 
summary(data.pca)

library(ggplot2)
library(factoextra)
library(ggpubr)
pdf(file="碎石图.PDF", height=5, width=6)
p<-fviz_screeplot(data.pca, addlabels = TRUE, main = NULL, 
                  ylim = c(0,80), barfill = "skyblue",
                  barcolor = "skyblue",ncp=10,linetype=2)+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 14), 
        plot.subtitle = element_text(size = 14), 
        plot.caption = element_text(size = 14), 
        strip.text = element_text(size = 14))
print(p)
dev.off()

library("corrplot") 
var = get_pca_var(data.pca)
pdf(file="contrib.PDF", height=5, width=6)
corrplot(var$contrib, is.corr=FALSE)
dev.off()

Type <- c(substr(rownames(data),11,17))
pdf(file="PCA_2d.PDF", height=5, width=6)
par(oma=c(0.5,0.5,0.5,0.5))
fviz_pca_biplot(data.pca, label = "var", habillage=Type,
                addEllipses=TRUE, ellipse.level=0.95,
                ggtheme = theme_minimal())+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12), 
        plot.subtitle = element_text(size = 12), 
        plot.caption = element_text(size = 12), 
        strip.text = element_text(size = 12))
dev.off()
