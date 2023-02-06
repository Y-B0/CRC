library(stringr)
library(clusterProfiler)
library(magrittr)
setwd("~/Desktop/project/TCGA/")
source("~/Desktop/project/TCGA/deg_script.R")

##TCGA
TCGA.READ.htseq_counts <- read.delim("~/Downloads/TCGA-READ.htseq_counts.tsv", header=T,row.names = 1)
samplename<-colnames(TCGA.READ.htseq_counts)
group<-data.frame(sample=samplename,
                  newname=ifelse(samplename%in%grep("01",grep("*\\.*\\.*\\.01",samplename,value = T),value = T),paste("T",samplename,sep = "_"),paste("C",samplename,sep = "_")),
                  group=ifelse(samplename%in%grep("01",grep("*\\.*\\.*\\.01",samplename,value = T),value = T),"T","C"))
group<-group[order(group$newname),]
TCGA.READ.htseq_counts<-TCGA.READ.htseq_counts[,group$sample]
##deg
deg_list<-limma_count2deg(TCGA.READ.htseq_counts,group$group,return.count = T)
deg<-deg_list[[1]]

deg$id<-str_split(deg$id,"\\.",simplify = T)[,1]

deg<-ensembl2symbol(deg)

volcano(deg,adj = "adj.P.Val")

save.image("READ.RData")

##heatmap
library(dplyr)
library(ComplexHeatmap)
logcpm<-deg_list[[2]]
et<-deg[abs(deg$log2FC)>=2&deg$Pvalue<0.01,] %>%.[order(.$log2FC),]
et<-rbind(head(et,100),tail(et,100))%>%.[order(.$log2FC,decreasing = T),]
rownames(logcpm)<-str_split(rownames(logcpm),"\\.",simplify = T)[,1]
et<-et[et$hgnc_symbol!="",]
logcpm<-logcpm[rownames(logcpm)%in%intersect(rownames(logcpm),et$id),]
rownames(logcpm)<-et$hgnc_symbol
#logcpm<-logcpm[str_split(rownames(logcpm),"\\.",simplify = T)[,1] %in% et$id,]%>%.[et$id,]

Heatmap(logcpm,row_dend_reorder = TRUE,row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),cluster_columns=T)
write.table(deg,file = "DEGs.xls",sep = "\t",quote = F,col.names = T)

##wgcna
deg<-deg[deg$adj.P.Val<0.05&abs(deg$log2FC)>0.585,]
deg<-deg[deg$hgnc_symbol!="",]
logcpm<-deg_list[[2]]
rownames(logcpm)<-str_split(rownames(logcpm),"\\.",simplify = T)[,1]
logcpm<-logcpm[rownames(logcpm)%in%deg$id,]
rownames(logcpm)<-deg$hgnc_symbol
feature<-as.data.frame.array(table(group$sample,group$group))


m.mad <- apply(logcpm,1,mad)
dataExpr <- logcpm[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]%>%t(.)
library(WGCNA)
powers<-1:50
nSamples=nrow(dataExpr)
ngene <- ncol(dataExpr)
type="signed"
corType = "bicor"
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
sft <- pickSoftThreshold(dataExpr, powerVector = powers)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "SFT, signed R^2", type = "n", 
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = 1:20, col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     type = "n", 
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, col = "red")
power<-sft[[1]]

sampleTree = hclust(dist(dataExpr), method = "average")

net = blockwiseModules(dataExpr, power = power, maxBlockSize = 10000,
                       TOMType = type,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       loadTOMs=TRUE,maxPOutliers=maxPOutliers,
                       saveTOMFileBase = paste0("result", ".tom"),
                       verbose = 3)
mergedColors = labels2colors(net$colors)
##plot cluster


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs
MEs_col = net$MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90,)

sampleName = rownames(dataExpr)
traitData = feature[match(sampleName, rownames(feature)), ]

if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               xLabelsAngle = 0,
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               setStdMargins = FALSE, textMatrix = textMatrix,
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### ???????????????????????????????????????

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# ???????????????????????????????????????

## ????????????????????????????????????????????????????????????????????????????????????????????????0-1?????????

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

##plot MM-GS
modNames = substring(names(MEs_col), 3)
setwd("./wgcna/GM_MM/")
for (i in modNames) {
  module = i
  for (j in colnames(traitData)) {
    pheno=j
    module_column = match(module, modNames)
    pheno_column = match(pheno,colnames(traitData))
    moduleGenes = moduleColors==module
    pdf(file = paste(module,pheno,".pdf",sep = "_"),width = 5,height = 5)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitCor[moduleGenes, pheno_column]),
                       xlab = paste("Module Membership in", module, "module and",pheno),
                       ylab = "Gene significance for body weight",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
  }
}
setwd("../../")


load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")
for (i in modNames) {
  module=i
  probes = colnames(dataExpr)
  inModule = (moduleColors==module)
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("./wgcna/network/",module,"model", ".edges.txt", sep=""),
                                 nodeFile = paste("./wgcna/network/",module,"model", ".nodes.txt", sep=""),
                                 weighted = TRUE, threshold = 0,
                                 nodeNames = modProbes, nodeAttr = moduleColors[inModule])
}
save.image("wgcna.RData")
##
deg<-deg[deg$adj.P.Val<0.01&abs(deg$log2FC)>2,]%>%.[.$hgnc_symbol!="",]

deg<-deg[order(deg$log2FC,decreasing = T),]

enrich_result<-enrich(deg$hgnc_symbol)

gene=bitr(deg$hgnc_symbol,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data_all <- merge(deg,gene,by.x="hgnc_symbol",by.y="SYMBOL")%>% arrange(desc(log2FC))
geneList = data_all$log2FC
names(geneList) <- data_all$ENTREZID

kegg_gmt<-read.gmt("~/Documents/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c2.cp.biocarta.v7.4.entrez.gmt")

gse.GO <- gseGO(
  geneList, 
  ont = "BP",  
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

gse.KEGG <- gseKEGG(geneList, 
                    organism = "hsa", 
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH")
library(enrichplot)

test<-gsea(deg)
write.csv(test[[1]]@result,"gse_GO.csv",row.names = F)
write.csv(test[[2]]@result,"gse_KEGG.csv",row.names = F)
gseaplot2(gse.KEGG, 1:5,ES_geom = "dot", base_size = 20, pvalue_table = T) 



###filter gene in model
load("wgcna.RData")
rm(list=ls()[!ls()%in%c("modTraitCor","geneModuleMembership","geneModuleMembershipA","geneTraitCor",
                        "geneTraitP","net")])
library(reshape2)
library(WGCNA)
gene<-as.data.frame(net$colors)
gene$`net$colors`<-paste("ME",gene$`net$colors`)
gene$color<-paste0("ME",labels2colors(as.numeric(str_replace_all(gene$`net$colors`,"ME",""))))
gene$gene<-names(net$colors)

module_g<-gene[,2:3]
module_g$color<-factor(module_g$color,level = names(table(module_g$color)[order(table(module_g$color))]))
module_g<-module_g[order(module_g$color),]
module_g<-split(module_g,module_g$color)
module_g<-lapply(module_g,function(x){x[,2]})
module_g<-do.call(cbind, lapply(lapply(module_g, unlist), `length<-`, max(lengths(module_g))))
module_g[is.na(module_g)]<-""
write.table(module_g,file = "module_gene.txt",sep = "\t",quote = F,row.names = F,col.names = T)

### enrich

enrich_list<-list()
cytoband_list<-list()
gene<-module_g
cytoband_gmt <- read.gmt("~/Documents/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c1.all.v7.4.entrez.gmt")
source("~/Documents/my_function/go_demo.R")
for (i in colnames(gene)) {
  key<-i
  genelist<-gene[,key]
  genelist<-bitr(genelist,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  try({
    dir.create(path=paste("enrich/",key,sep = ""),recursive = T)
    setwd(paste("enrich/",key,sep = ""))
    tmp<-go_demo(genesymbol = genelist$SYMBOL,plot = F)
    enrich_list[[key]]<-tmp
    cytoband<-enricher(genelist$ENTREZID,TERM2GENE = cytoband_gmt)
    cytoband_list[[key]]<-cytoband@result
    write.table(as.data.frame(cytoband@result),file = "cytoband.txt",row.names = F,col.names = T,quote = F,sep = "\t")
  })
  setwd("../../")
}
save(enrich_list,file="enrich/enrich_model.RData")
save(cytoband_list,file = "enrich/cytoband_module.RData")
##process enrich result
load("enrich/enrich_model.RData")
for (i in 1:length(enrich_list)) {
  try({
    enrich_list[[i]]$method<-names(enrich_list[i])
    enrich_list[[i]]$model<-names(enrich_list[i])
    #enrich_list[[i]]<-enrich_list[[i]][enrich_list[[i]]$ONTOLOGY=="BP",]
    enrich_list[[i]]<-enrich_list[[i]][enrich_list[[i]]$pvalue<0.05&enrich_list[[i]]$qvalue<0.05,]
    enrich_list[[i]]<-split(enrich_list[[i]],enrich_list[[i]]$ONTOLOGY)%>%lapply(.,head,n=10)%>%Reduce(rbind,.)
    #enrich_list[[i]]<-head(enrich_list[[i]],n=20)
  })
}
enrich_list<-Reduce(rbind,enrich_list)
##go 9:10
enrich_list<-enrich_list[enrich_list$model=="MEbrown",]
enrich_list<-enrich_list[order(enrich_list$Count,decreasing = T),]
enrich_list<-enrich_list[order(enrich_list$ONTOLOGY,decreasing = F),]
enrich_list$Description <- factor(enrich_list$Description,levels = rev(enrich_list$Description))

ggplot(data = enrich_list, # ?????????????????????
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # ?????????????????????????????????
  geom_bar(stat = "identity",width = 0.9)+ # ??????????????????????????????
  coord_flip()+theme_bw()+ # ????????????????????????????????????
  scale_x_discrete(labels = function(x) str_wrap(x,width = 80))+ # ??????term?????????????????????
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # ??????????????????????????????
  theme(axis.title = element_text(size = 13), # ?????????????????????
        axis.text = element_text(size = 11), # ?????????????????????
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # ????????????
        legend.title = element_text(size = 13), # ??????????????????
        legend.text = element_text(size = 11), # ??????????????????
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))+ # ?????????
  scale_fill_aaas()

## TF

TFs <- read.delim("Integrated_meanRank.tsv")%>%head(.,n=10)
tf_list<-separate_rows(TFs, Overlapping_Genes,sep = ",")%>%as.data.frame(.)
tmp<-as.data.frame(table(tf_list$TF))%>%.[.$Freq!=0,]
tmp<-tmp[order(tmp$Freq,decreasing = T),]
tf_list$TF<-factor(tf_list$TF,levels = tmp$Var1)
tf_list<-tf_list[order(tf_list$TF),]
tf_list<-split(tf_list,tf_list$TF)
tf_list<-lapply(tf_list,function(x){x[,6]})
hubgene<-Reduce(intersect,tf_list)%>%as.data.frame(.)
tf_list<-do.call(cbind, lapply(lapply(tf_list, unlist), `length<-`, max(lengths(tf_list))))
tf_list[is.na(tf_list)]<-""
write.table(tf_list,file = "tf_gene.txt",row.names = F,col.names = T,sep = "\t",quote = F)
write.table(hubgene,file = "hubgene.txt",row.names = F,col.names = T,sep = "\t",quote = F)
tmp$proportion<-tmp$Freq/1082;tmp$proportion<-round(tmp$proportion*100,0)
tmp$Var1<-as.character(tmp$Var1)
tmp$proportion<-tmp$proportion*0.01
tmp$Var1<-factor(tmp$Var1,levels = tmp$Var1)
ggplot(tmp)+geom_bar(aes(x=Var1,y=proportion),position = "dodge",width = 0.6,stat = "identity")+theme_bw()+ylab("Rate")+xlab("")+
  theme(text=element_text(size=20),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+
  geom_text(aes(x = Var1, y= 0.5,label = Freq),size=4.5,check_overlap = T)


##
library(psych)
library(ggcorrplot)
load("READ.RData")
rm(list = ls()[!ls()%in%c("logcpm","group","deg")])
rownames(logcpm)<-str_split(rownames(logcpm),"\\.",simplify = T)[,1]
deg<-deg[deg$hgnc_symbol!="",]
deg<-deg[!duplicated(deg$hgnc_symbol),]
deg<-deg[!duplicated(deg$id),]
deg<-deg[deg$id%in%intersect(rownames(logcpm),deg$id),]
logcpm<-logcpm[rownames(logcpm)%in%intersect(rownames(logcpm),deg$id),]
rownames(logcpm)<-deg$hgnc_symbol

TFs <- read.delim("Integrated_meanRank.tsv")%>%head(.,n=10)
ICD_gene <- read_excel("ICD_gene.xlsx")
hubgene<-read.table("hubgene.txt",header = T) 
cor<-corr.test(t(logcpm[intersect(rownames(logcpm),ICD_gene$`ICD (10.3389/fgene.2022.1069921)`),]),
               t(logcpm[c(intersect(rownames(logcpm),TFs$TF),intersect(rownames(logcpm),hubgene$hubgene)),]),method = "spearman")

ggcorrplot(cor$r,p.mat = cor$p,hc.order = F,outline.color = "white",
              ggtheme = theme_bw(),sig.level = 0.05,insig = "blank",lab_size = 3,digits = 1,
              colors = c("#0081A7","#FDFCDC","#F07167"),lab = T)+
  theme(axis.text.x = element_text(color = "red"),axis.text.y = element_text(color = rep(c("#C24294","#4EADEA"),c(10,35))),axis.text.x.top = element_text(hjust = 0,vjust = 0.5,angle = 90))+
  scale_x_discrete(position = "top")

string_gene_list<-data.frame(gene=c(rownames(cor$r),colnames(cor$r)),type=rep(c("ICD","hubgene","TF"),times=c(54,35,10)))
write.table(string_gene_list,file = "string_gene_list.txt",row.names = F,col.names = T,sep = "\t",quote = F)



## survival

library(survival)
library(survminer)
load("READ.RData")
rm(list = ls()[!ls()%in%c("logcpm","group","deg")])
rownames(logcpm)<-str_split(rownames(logcpm),"\\.",simplify = T)[,1]
deg<-deg[deg$hgnc_symbol!="",]
deg<-deg[!duplicated(deg$hgnc_symbol),]
deg<-deg[!duplicated(deg$id),]
deg<-deg[deg$id%in%intersect(rownames(logcpm),deg$id),]
logcpm<-logcpm[rownames(logcpm)%in%intersect(rownames(logcpm),deg$id),]
rownames(logcpm)<-deg$hgnc_symbol
TFs <- read.delim("Integrated_meanRank.tsv")%>%head(.,n=10)

logcpm<-logcpm[intersect(TFs$TF,rownames(logcpm)),group$sample[group$group=="T"]]
logcpm<-t(logcpm)

## lasso module

##EaSleR




