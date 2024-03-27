# Script to plot using heatmaps the proteins of GradR in Nostoc

data<-read.table("Full_dataset_V3.csv",header=TRUE,sep=";")
nrow(data)

# Individual tables for each gradient
table.rnase1<-data[,seq(11,113,by=6)]*100
rownames(table.rnase1)<-data$Gene
table.rnase1["patU3",]
table.rnase2<-data[,seq(12,114,by=6)]*100
rownames(table.rnase2)<-data$Gene
table.rnase3<-data[,seq(13,115,by=6)]*100
rownames(table.rnase3)<-data$Gene
table.ctrl1<-data[,seq(14,116,by=6)]*100
rownames(table.ctrl1)<-data$Gene
table.ctrl2<-data[,seq(15,117,by=6)]*100
rownames(table.ctrl2)<-data$Gene
table.ctrl3<-data[,seq(16,118,by=6)]*100
rownames(table.ctrl3)<-data$Gene

# Filter the protein. Only keep proteins that are expressed in the proportions in at least one ctrl and one rnase

proteins.identified<-which((rowSums(table.ctrl1) + rowSums(table.ctrl2) + rowSums(table.ctrl3)) > 0 & (rowSums(table.rnase1) + rowSums(table.rnase2) + rowSums(table.rnase3)) > 0)
length(proteins.identified)

length(proteins.identified)# 2638 proteins identified just looking to identification but only 1848 with expression

protein.identification.table<-data.frame(ctrl1=rep(0,length(proteins.identified)),ctrl2=rep(0,length(proteins.identified)),ctrl3=rep(0,length(proteins.identified))
                                         ,rnase1=rep(0,length(proteins.identified)),rnase2=rep(0,length(proteins.identified),rnase3=rep(0,length(proteins.identified))))
rownames(protein.identification.table)<-names(proteins.identified)

for(i in 1:nrow(protein.identification.table))
{
  if(rowSums(table.ctrl1)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"ctrl1"]<-1
  }
  if(rowSums(table.ctrl2)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"ctrl2"]<-1
  }
  if(rowSums(table.ctrl3)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"ctrl3"]<-1
  }
  if(rowSums(table.rnase1)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"rnase1"]<-1
  }
  if(rowSums(table.rnase2)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"rnase2"]<-1
  }
  if(rowSums(table.rnase3)[names(proteins.identified[i])] > 0)
  {
    protein.identification.table[names(proteins.identified[i]),"rnase3"]<-1
  }
}


protein.identification.table <- rapply(protein.identification.table, f=function(x) ifelse(is.na(x),0,x), how="replace") 

protein.identification.table[["n_ctrl"]]<-rep(0,nrow(protein.identification.table))
protein.identification.table[["n_rnase"]]<-rep(0,nrow(protein.identification.table))

for(i in 1:nrow(protein.identification.table))
{
  protein.identification.table$n_ctrl[i]<-sum(protein.identification.table[i,1:3])
  protein.identification.table$n_rnase[i]<-sum(protein.identification.table[i,4:6])
}

# Create the combined dataset for our proteins

# Create a dataframe of expression with the proteins used for my analysis
table.ctrl1x<-table.ctrl1
nrow(table.ctrl1x)
rownames(table.ctrl1x)<-paste(rownames(table.ctrl1x),"ctrl1",sep="_")
table.ctrl1x$Gene<-data$Gene
colnames(table.ctrl1x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.ctrl1x<-table.ctrl1x[proteins.identified,]

table.ctrl2x<-table.ctrl2
rownames(table.ctrl2x)<-paste(rownames(table.ctrl2x),"ctrl2",sep="_")
table.ctrl2x$Gene<-data$Gene
colnames(table.ctrl2x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.ctrl2x<-table.ctrl2x[proteins.identified,]

table.ctrl3x<-table.ctrl3
rownames(table.ctrl3x)<-paste(rownames(table.ctrl3x),"ctrl3",sep="_")
table.ctrl3x$Gene<-data$Gene
colnames(table.ctrl3x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.ctrl3x<-table.ctrl3x[proteins.identified,]

table.rnase1x<-table.rnase1
rownames(table.rnase1x)<-paste(rownames(table.rnase1x),"rnase1",sep="_")
table.rnase1x$Gene<-data$Gene
colnames(table.rnase1x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.rnase1x<-table.rnase1x[proteins.identified,]

table.rnase2x<-table.rnase2
rownames(table.rnase2x)<-paste(rownames(table.rnase2x),"rnase2",sep="_")
table.rnase2x$Gene<-data$Gene
colnames(table.rnase2x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.rnase2x<-table.rnase2x[proteins.identified,]

table.rnase3x<-table.rnase3
rownames(table.rnase3x)<-paste(rownames(table.rnase3x),"rnase3",sep="_")
table.rnase3x$Gene<-data$Gene
colnames(table.rnase3x)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
table.rnase3x<-table.rnase3x[proteins.identified,]

table.expression.combined.normalized<-rbind(table.ctrl1x,table.ctrl2x,table.ctrl3x,table.rnase1x,table.rnase2x,table.rnase3x)
colnames(table.expression.combined.normalized)<-c("fraction1_2",paste("fraction",3:19,sep=""),"Gene")
head(table.expression.combined.normalized)

write.table(table.expression.combined.normalized,"table_expression_combined_proportions.csv",sep="\t")

# Check the correlation between samples

matrix.correlations<-matrix(rep(NA,nrow(table.ctrl1x))*9,nrow=nrow(table.ctrl1x),ncol=9)
rownames(matrix.correlations)<-table.ctrl1x$Gene
colnames(matrix.correlations)<-c("ctrl1_vs_ctrl2","ctrl1_vs_ctrl3","ctrl2_vs_ctrl3","rnase1_vs_rnase2","rnase1_vs_rnase3","rnase2_vs_rnase3","ctrl1_vs_rnase1","ctrl2_vs_rnase2","ctrl3_vs_rnase3")
head(matrix.correlations)

for(i in 1:nrow(table.ctrl1x))
{
  matrix.correlations[i,1]<-cor(t(table.ctrl1x[i,1:18]),t(table.ctrl2x[i,1:18]),method="spearman")
  matrix.correlations[i,2]<-cor(t(table.ctrl1x[i,1:18]),t(table.ctrl3x[i,1:18]),method="spearman")
  matrix.correlations[i,3]<-cor(t(table.ctrl2x[i,1:18]),t(table.ctrl3x[i,1:18]),method="spearman")
  matrix.correlations[i,4]<-cor(t(table.rnase1x[i,1:18]),t(table.rnase2x[i,1:18]),method="spearman")
  matrix.correlations[i,5]<-cor(t(table.rnase1x[i,1:18]),t(table.rnase3x[i,1:18]),method="spearman")
  matrix.correlations[i,6]<-cor(t(table.rnase2x[i,1:18]),t(table.rnase3x[i,1:18]),method="spearman")
  matrix.correlations[i,7]<-cor(t(table.ctrl1x[i,1:18]),t(table.rnase1x[i,1:18]),method="spearman")
  matrix.correlations[i,8]<-cor(t(table.ctrl2x[i,1:18]),t(table.rnase2x[i,1:18]),method="spearman")
  matrix.correlations[i,9]<-cor(t(table.ctrl3x[i,1:18]),t(table.rnase3x[i,1:18]),method="spearman")
  
}

write.table(matrix.correlations,"matrix_correlations.csv",sep="\t")

data.frame.correlations<-as.data.frame(matrix.correlations)
head(data.frame.correlations)

library(ggplot2)

ggplot(data.frame.correlations, aes(x=ctrl1_vs_ctrl2)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=ctrl1_vs_ctrl3)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=ctrl2_vs_ctrl3)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=rnase1_vs_rnase2)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=rnase1_vs_rnase3)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=rnase2_vs_rnase3)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=ctrl1_vs_rnase1)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=ctrl2_vs_rnase2)) + geom_histogram(binwidth=0.01,color="black",fill="white")
ggplot(data.frame.correlations, aes(x=ctrl3_vs_rnase3)) + geom_histogram(binwidth=0.01,color="black",fill="white")



mean(matrix.correlations[,1],na.rm=TRUE)
mean(matrix.correlations[,2],na.rm=TRUE)
mean(matrix.correlations[,3],na.rm=TRUE)
mean(matrix.correlations[,4],na.rm=TRUE)
mean(matrix.correlations[,5],na.rm=TRUE)
mean(matrix.correlations[,6],na.rm=TRUE)
mean(matrix.correlations[,7],na.rm=TRUE)
mean(matrix.correlations[,8],na.rm=TRUE)
mean(matrix.correlations[,9],na.rm=TRUE)

rm(table.data.ctrl1.export,table.data.ctrl2.export,table.data.ctrl3.export,table.data.rnase1.export,table.data.rnase2.export,table.data.rnase3.export)


#Creation of a dataframe for the mean of several replicates

table.mean.expression.ctrl<-table.ctrl1x
rownames(table.mean.expression.ctrl)<-paste(table.mean.expression.ctrl$Gene,"ctrl",sep="_")
table.mean.expression.ctrl[,1:18]<-0

table.mean.expression.rnase<-table.rnase1x
rownames(table.mean.expression.rnase)<-paste(table.mean.expression.rnase$Gene,"rnase",sep="_")
table.mean.expression.rnase[,1:18]<-0

for(i in 1:nrow(table.mean.expression.ctrl))
{
    table.mean.expression.ctrl[i,1:18]<-(table.ctrl1x[i,1:18] + table.ctrl2x[i,1:18] + table.ctrl3x[i,1:18])/3
    table.mean.expression.rnase[i,1:18]<-(table.rnase1x[i,1:18] + table.rnase2x[i,1:18] + table.rnase3x[i,1:18])/3
}

matrix.mean.expression.ctrl.normalized<-as.matrix(table.mean.expression.ctrl[,1:18]*100 / rowSums(table.mean.expression.ctrl[,1:18]))
rownames(matrix.mean.expression.ctrl.normalized)<-table.mean.expression.ctrl$Gene
matrix.mean.expression.rnase.normalized<-as.matrix(table.mean.expression.rnase[,1:18]*100 / rowSums(table.mean.expression.rnase[,1:18]))
rownames(matrix.mean.expression.rnase.normalized)<-table.mean.expression.rnase$Gene

for(i in 1:nrow(table.mean.expression.ctrl))
{
  table.mean.expression.ctrl[i,1:18]<-matrix.mean.expression.ctrl.normalized[i,]
  table.mean.expression.rnase[i,1:18]<-matrix.mean.expression.rnase.normalized[i,]
}


table.mean.expression.ctrl$Found_in_replicates<-rep("NA",nrow(table.mean.expression.ctrl))
table.mean.expression.rnase$Found_in_replicates<-rep("NA",nrow(table.mean.expression.rnase))

for(i in 1:nrow(table.mean.expression.ctrl))
{
  table.mean.expression.ctrl$Found_in_replicates[i]<-sum(protein.identification.table[table.mean.expression.ctrl$Gene[i],c("ctrl1","ctrl2","ctrl3")])
  table.mean.expression.rnase$Found_in_replicates[i]<-sum(protein.identification.table[table.mean.expression.ctrl$Gene[i],c("rnase1","rnase2","rnase3")])
  
}

table.mean.expression.combined<-rbind(table.mean.expression.ctrl,table.mean.expression.rnase)
write.table(table.mean.expression.combined,"table_mean_expression_combined.csv",sep="\t")

# Load the limma data

limma_data<-read.table("limma_results_proportions_V3.csv",header=TRUE,sep=",")
limma_data<-limma_data[,c("Gene","Protein.ID","Protein.Description","sample","logFC","fdr.limma","qval.fdrtool","hit_annotation","hit_annotation_method")]

matrix.foldchange<-matrix(rep(0,length(proteins.identified)*18),nrow=length(proteins.identified),ncol=18)
rownames(matrix.foldchange)<-names(proteins.identified)
colnames(matrix.foldchange)<-c("fraction01_2",paste("fraction0",3:9,sep=""),paste("fraction",10:19,sep=""))
head(matrix.foldchange)

matrix.pvalue<-matrix(rep(0,length(proteins.identified)*18),nrow=length(proteins.identified),ncol=18)
rownames(matrix.pvalue)<-names(proteins.identified)
colnames(matrix.pvalue)<-c("fraction01_2",paste("fraction0",3:9,sep=""),paste("fraction",10:19,sep=""))
head(matrix.pvalue)

names(proteins.identified)

for(i in 1:nrow(limma_data))
{
  print(i)
  if(sum(names(proteins.identified) == limma_data$Gene[i]) > 0 )
  {
    matrix.foldchange[limma_data$Gene[i],limma_data$sample[i]]<-limma_data$logFC[i]
    matrix.pvalue[limma_data$Gene[i],limma_data$sample[i]]<-limma_data$qval.fdrtool[i]
  }
}

table.foldchange<-as.data.frame(matrix.foldchange)
table.foldchange$Gene<-rownames(table.foldchange)
table.foldchange$Protein.ID<-NA
table.foldchange$Protein.Description<-NA
table.foldchange$n_ctrl<-NA
table.foldchange$n_rnase<-NA
table.foldchange$n_hits<-0
table.foldchange$hits<-NA
table.foldchange$n_candidates<-0
table.foldchange$candidates<-NA

for(i in 1:nrow(table.foldchange))
{
  table.foldchange$n_ctrl[i]<-protein.identification.table[rownames(table.foldchange)[i],"n_ctrl"]
  table.foldchange$n_rnase[i]<-protein.identification.table[rownames(table.foldchange)[i],"n_rnase"]
}

for(i in 1:nrow(limma_data))
{
  if(sum(rownames(table.foldchange) == limma_data$Gene[i]) > 0 )
  {
    table.foldchange[limma_data$Gene[i],"Protein.ID"]<-limma_data$Protein.ID[i]
    table.foldchange[limma_data$Gene[i],"Protein.Description"]<-limma_data$Protein.Description[i]
    table.foldchange[limma_data$Gene[i],"n_hits"]<-sum(limma_data[limma_data$Gene[i] == limma_data$Gene,"hit_annotation"] == "hit")
    table.foldchange[limma_data$Gene[i],"n_candidates"]<-sum(limma_data[limma_data$Gene[i] == limma_data$Gene,"hit_annotation"] == "candidate")
    limma_my_protein<-limma_data[limma_data$Gene[i] == limma_data$Gene,]
    if(sum(limma_my_protein$hit_annotation == "hit") > 0 | sum(limma_my_protein$hit_annotation == "candidate") > 0)
    {
      if(sum(limma_my_protein$hit_annotation == "hit") > 0)
      {
        table.foldchange[limma_data$Gene[i],"hits"]<-paste(limma_my_protein[limma_my_protein$hit_annotation == "hit","sample"],collapse=" ")
      }
      if(sum(limma_my_protein$hit_annotation == "candidate") > 0)
      {
        table.foldchange[limma_data$Gene[i],"candidates"]<-paste(limma_my_protein[limma_my_protein$hit_annotation == "candidate","sample"],collapse=" ")
      } 
    }
  }
}

colnames(table.foldchange)[1:18]<-c("FC_fraction01_2",paste("FC_fraction0",3:9,sep=""),paste("FC_fraction",10:19,sep=""))

table.foldchange[["n_hits_plus_candidates"]]<- table.foldchange$n_hits + table.foldchange$n_candidates

colnames(matrix.pvalue)<-c("p_fraction01_2",paste("p_fraction0",3:9,sep=""),paste("p_fraction",10:19,sep=""))

table.final<-cbind(table.foldchange,as.data.frame(matrix.pvalue))

write.table(table.final,"table.final.csv",sep="\t")

# Load the SVM information

table.SVM<-read.table("SVM_prediction_Halie.csv",sep="\t",header=TRUE)
nrow(table.SVM)
rownames(table.SVM)<-table.SVM$Symbol

table.final$SVM<-table.SVM[rownames(table.final),"Score"]

write.table(table.final,"table.final2.csv",sep="\t")

# GO enrichment of the candidates

table.final<-read.table("table.final.csv",header=TRUE,sep="\t",row.names=1)
head(table.final)
nrow(table.final)

candidates<-rownames(table.final)[table.final$n_hits_plus_candidates > 0]
length(candidates)

library(clusterProfiler)
library(pathview)
library(enrichplot)
library(GOplot)
library(DOSE)
library(topGO)
library(org.N7120.eg.db)


all_symbol<-rownames(table.final)

ego<-enrichGO(gene= candidates, OrgDb = org.N7120.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff= 0.05, universe = all_symbol, keyType = "SYMBOL")

ego.res<-as.data.frame(ego)
head(ego.res)

write.table(ego.res,"GO_enrichment_candidates.csv",sep="\t")


# InterPro Enrichment


significant_interpro<-function(candidates,background,type)
{
  Interpro_ID012<-AnnotationDbi::select(org.N7120.eg.db, keys=candidates, columns=c("InterPro","Interpro_type","Interpro_name"), keytype="SYMBOL")
  
  Interpro_no_NAs<-Interpro_ID012[!is.na(Interpro_ID012$Interpro_type),]
  Interpro_only_type<-Interpro_no_NAs[Interpro_no_NAs$Interpro_type == type,]
  
  unique_only_type<-unique(Interpro_only_type$InterPro)
  type_p_values<-data.frame(InterProID="InterProID",P.value="p.value",Proportion.candidates="Proportion_candidates",Proportion.background="Proportion_background",Gene="Symbol")
  
  Interpro_background<-AnnotationDbi::select(org.N7120.eg.db, keys=background, columns=c("InterPro","Interpro_type","Interpro_name"), keytype="SYMBOL")
  Interpro_background_no_NAs<-Interpro_background[!is.na(Interpro_background$Interpro_type),]
  Interpro_background_only_type<-Interpro_background_no_NAs[Interpro_background_no_NAs$Interpro_type == type,]
  
  
  for (i in 1:length(unique_only_type))
  {
    print(i)
    q=sum(Interpro_only_type$InterPro == unique_only_type[i])
    m=sum(Interpro_background_only_type$InterPro == unique_only_type[i])
    n=length(Interpro_background_only_type$InterPro)- sum(Interpro_background_only_type$InterPro == unique_only_type[i])
    k=length(Interpro_only_type$InterPro)
    p.value<-phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
    proportion_candidates<- paste(as.character(q),"/", as.character(k),sep="")
    proportion_background<- paste(as.character(m),"/", as.character(length(Interpro_background_only_type$InterPro)),sep="")
    test<-data.frame(InterProID=unique_only_type[i],P.value=p.value,Proportion.candidates=proportion_candidates,Proportion.background=proportion_background,Gene="")
    type_p_values<-rbind(type_p_values,test)
  }
  
  type_p_values<-type_p_values[2:nrow(type_p_values),]
  type_p_values$Name<-rep("NA",nrow(type_p_values))
  type_p_values$P.value.adjust<-p.adjust(type_p_values$P.value,"BH")
  
  for(i in 1:nrow(type_p_values))
  {
    for(j in 1:nrow(Interpro_no_NAs))
    {
      if(type_p_values$InterPro[i] == Interpro_no_NAs$InterPro[j])
      {
        type_p_values$Name[i]<-Interpro_no_NAs$Interpro_name[j]
        type_p_values$Gene[i]<-paste(type_p_values$Gene[i],"/",Interpro_no_NAs$SYMBOL[j],sep="")
      }
    }
  }
  
  return(type_p_values)
}

enrichment_interpro<-significant_interpro(candidates,all_symbol,"Homologous_superfamily")
enrichment_interpro

# Nothing enriched

## Clustering analysis
# Load a table with the average of replicates only for the control

table_clustering<-read.table("table_clustering_ctrl.csv",header=TRUE,sep="\t",row.names=1) # It is the table.mean.expression.combined but only the ctrl samples

################################
################# Mfuzz approach

## Clustering with Mfuzz
library(Mfuzz)
library(limma)
library(factoextra)

eset <- new("ExpressionSet",exprs=as.matrix(table_clustering))

cyano.s <- Mfuzz::standardise(eset)

m1<-mestimate(cyano.s)
m1

# Calculate the best cluster number from 2 to 30

vector_silhouette2<-c()

for(i in 2:30)
{
  print(i)
  cluster.number<-i
  set.seed(2)
  cl <- mfuzz(eset,c=cluster.number,m=1.14)
  acore_list<-acore(eset,cl=cl,min.acore=0.85)
  cluster.groups<- data.frame(Gene = rownames(table_clustering),
                              soft.cluster.group = rep(NA,nrow(table_clustering)))
  for(j in 1:nrow(cluster.groups))
  {
    for(z in 1:cluster.number)
    {
      if(sum(cluster.groups$Gene[j] == acore_list[[z]]$NAME) > 0)
      {
        cluster.groups[j,"soft.cluster.group"]<-z
      }
    }
  }
  factor_soft_cluster<-cluster.groups$soft.cluster.group
  names(factor_soft_cluster)<-cluster.groups$Gene
  factor_soft_cluster<-factor_soft_cluster[!is.na(factor_soft_cluster)]
  soft.sil<-summary(cluster::silhouette(factor_soft_cluster, dist(table_clustering[names(factor_soft_cluster),])))[["avg.width"]]
  vector_silhouette2<-c(vector_silhouette2,soft.sil)
}

plot(2:30,vector_silhouette2,type="o",col="blue",ylim=c(0.2,0.5),pch=0,xlab="Number of clusters",ylab="Silhouette")


# 19 clusters

cluster.number<-19

set.seed(2)
cl <- mfuzz(eset,c=cluster.number,m=1.14) #I have selected 19 because was good number for kmeans

mfuzz.plot2(cyano.s,cl=cl,min.mem= 0.85,mfrow=c(4,5)) #Command to visualize clusters

acore_list<-acore(eset,cl=cl,min.acore=0.9)

# To facilitate the interpretation of the reader I have reordered the clusters, in this order: 7-2-6-19-9-10-17--5-13-18-4-11-12-3-15-8-1-14-16

new_acore_list<-acore_list[c(7,2,6,19,9,10,17,5,13,18,4,11,12,3,15,8,1,14,16)]

cluster.groups<- data.frame(Gene = rownames(table_clustering),
                            soft.cluster.group = rep(NA,nrow(table_clustering)),confidence=rep(NA,nrow(table_clustering)))
cluster

for(i in 1:nrow(cluster.groups))
{
  for(j in 1:cluster.number)
  {
    if(sum(cluster.groups$Gene[i] == new_acore_list[[j]]$NAME) > 0)
    {
      cluster.groups[i,"soft.cluster.group"]<-j
      cluster.groups[i,"confidence"]<-new_acore_list[[j]][cluster.groups$Gene[i],"MEM.SHIP"]
    }
  }
}

write.table(cluster.groups,"mfuzz_19_09confidence_reordered.csv",sep="\t")



# Calculate silhouette for my soft clustering

factor_soft_cluster<-cluster.groups$soft.cluster.group
names(factor_soft_cluster)<-cluster.groups$Gene

factor_soft_confidence<-cluster.groups$confidence
names(factor_soft_confidence)<-cluster.groups$Gene

factor_soft_cluster<-factor_soft_cluster[!is.na(factor_soft_cluster)]
factor_soft_confidence<-factor_soft_confidence[!is.na(factor_soft_confidence)]
length(factor_soft_cluster)

soft.sil<-cluster::silhouette(factor_soft_cluster, dist(table_clustering[names(factor_soft_cluster),]))
summary(soft.sil)[["avg.width"]]
fviz_silhouette(soft.sil)

## PCA of all proteins with the ctrl samples

clust.data_PCA <- prcomp(table_clustering, scale = FALSE)
perc_var <-
  round(100 * clust.data_PCA$sdev ^ 2 /
          sum(clust.data_PCA$sdev ^ 2), 1)
PCA_clust.data <-
  data.frame(PC1 = clust.data_PCA$x[, 1],
             PC2 = clust.data_PCA$x[, 2],
             PC3 = clust.data_PCA$x[, 3])
ggplot(data = PCA_clust.data, aes(PC1, PC2, colour = PC3)) +
  geom_point() +
  theme_bw(base_size = 12) +
  scale_colour_gradientn(colours = c("black", "#377eb8", "#984ea3", 
                                     "#e41a1c", "#ff7f00", "#ffff33")) +
  ggtitle("PCA clustering data",
          subtitle = paste("PC3 - explaining", perc_var[3], "% of variability")) +
  xlab(paste("PC1 - explaining", perc_var[1], "% of variability")) +
  ylab(paste("PC2 - explaining", perc_var[2], "% of variability"))


# Plot soft cluster information in PCA
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)

cluster.groups<-cluster.groups[,c(1,2)]

m_cluster.groups <- cluster.groups %>%
  group_by(Gene) %>%
  gather(-Gene, key = "cluster.group", value = "cluster") %>%
  mutate(cluster.group = gsub(".cluster.group", "", cluster.group))

PCA_clust.data$Gene <- rownames(PCA_clust.data)
PCA_clust.data <- left_join(PCA_clust.data, m_cluster.groups)
PCA_clust.data$cluster <- factor(PCA_clust.data$cluster)
ggplot(data = PCA_clust.data, aes(PC1, PC2, colour = cluster)) +
  geom_point() +
  theme_bw(base_size = 12) +
  facet_wrap(~ cluster.group) +
  geom_text(aes(label = Gene), size = 0.2, colour = "black")

# Adjust table of clustering for plotting

table_clustering<-table_clustering[names(factor_soft_cluster),]#Only clustered data
nrow(table_clustering)
table_clustering$soft.cluster.group<-factor_soft_cluster
table_clustering$confidence<-factor_soft_confidence
table_clustering$Gene<-rownames(table_clustering)
m_clust.data<-reshape2::melt(table_clustering,id.vars=c("Gene","soft.cluster.group","confidence"))
head(m_clust.data)
colnames(m_clust.data)<-c("Gene","soft.cluster.group","confidence","sample","mean.value")
head(m_clust.data)

customPlot <- list(
  theme_bw(base_size = 12), 
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8)), 
  scale_fill_brewer(palette = "Set1"), 
  scale_colour_brewer(palette = "Set1")
)


gr.width <- 4 + ceiling(sqrt(cluster.number)) * 3
gr.height <- 2.5 + floor(sqrt(cluster.number)) * 3

ggplot(data = m_clust.data, aes(sample, mean.value)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = paste(Gene, soft.cluster.group)), alpha = 0.1) +
  geom_smooth(fun.data = "mean_se", stat = "summary", aes(group = soft.cluster.group)) +
  facet_wrap(~ soft.cluster.group) +
  customPlot +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("soft clustering") +
  ylab("% protein")
ggsave(file.path(paste0("Clustering_line_plot_soft_", cluster.number, "_cluster_", nrow(table_clustering), "_proteins_", ".pdf")), 
       width = gr.width, height = gr.height)

# Cluster the log2 foldchange data
table.final<-read.table("table.final.csv",sep="\t",header=TRUE,row.names=1)
head(table.final)

shifting_proteins<-read.table("333_shifting_proteins.txt",header=TRUE)$ID
table.clustering.log2<-table.final[shifting_proteins,1:18]
head(table.clustering.log2)
table.clustering.log2[is.na(table.clustering.log2)]<-0

#PCA of the log2 samples

clust.data_PCA <- prcomp(table.clustering.log2, scale = FALSE)
perc_var <-
  round(100 * clust.data_PCA$sdev ^ 2 /
          sum(clust.data_PCA$sdev ^ 2), 1)
PCA_clust.data <-
  data.frame(PC1 = clust.data_PCA$x[, 1],
             PC2 = clust.data_PCA$x[, 2],
             PC3 = clust.data_PCA$x[, 3])
ggplot(data = PCA_clust.data, aes(PC1, PC2, colour = PC3)) +
  geom_point() +
  theme_bw(base_size = 12) +
  scale_colour_gradientn(colours = c("black", "#377eb8", "#984ea3", 
                                     "#e41a1c", "#ff7f00", "#ffff33")) +
  ggtitle("PCA clustering data",
          subtitle = paste("PC3 - explaining", perc_var[3], "% of variability")) +
  xlab(paste("PC1 - explaining", perc_var[1], "% of variability")) +
  ylab(paste("PC2 - explaining", perc_var[2], "% of variability"))

# Check the best cluster number and clustered

library(NbClust)
library(factoextra)
library(gplots)
fviz_nbclust(table.clustering.log2, FUNcluster=hcut, method="silhouette")+theme_classic()

cluster.number<-2 # According to the silhouette analysis k=2 is the best option.

mypalette <- colorRampPalette(c("#006CD1","white","#994F00"))(n = 100)
index<-heatmap.2(as.matrix(table.clustering.log2),Colv="FALSE",dendrogram="row",scale="none",density.info="none",trace="none",key=TRUE, col= mypalette,cexRow = 0.1,cexCol = 0.8,na.color="grey")

hclust.fit<-as.hclust(index$rowDendrogram)

cluster.groups_hclust <- data.frame(Gene = names(cutree(hclust.fit, k=cluster.number)),
             hclust.cluster.group = cutree(hclust.fit, k = cluster.number))

write.table(cluster.groups_hclust,"hclust_2_clusters_log2.csv",sep="\t")

# Heatmap of protein distributions

library(gplots)

heatmap.plot.protein<-function(protein)
{
  my_protein<-protein
  matrix_heatmap<-as.matrix(table.expression.combined.normalized[table.expression.combined.normalized$Gene == my_protein,c(1:18)])
  names<-c()
  
  if(length(strsplit(rownames(matrix_heatmap)[1],split="_")[[1]]) == 2) # Chequear esto
  {
    for(i in 1:nrow(matrix_heatmap))
    {
      names<-c(names,strsplit(rownames(matrix_heatmap)[i],split="_")[[1]][2])
    }
  }
  
  if(length(strsplit(rownames(matrix_heatmap)[1],split="_")[[1]]) == 3) # Chequear esto
  {
    for(i in 1:nrow(matrix_heatmap))
    {
      names<-c(names,strsplit(rownames(matrix_heatmap)[i],split="_")[[1]][3])
    }
  }  
  
  rownames(matrix_heatmap)<-names
  mypalette <- colorRampPalette(c("white","black"))(n = 100)
  heatmap.2(matrix_heatmap,main=my_protein,Rowv="FALSE",Colv="FALSE",dendrogram="none",scale="none",density.info="none",trace="none",key=TRUE, col= mypalette,cexRow = 0.8)
  
}

heatmap.plot.multiple.proteins<-function(proteins,treatment)
{
  
  if(treatment == "ctrl")
  {
    matrix_heatmap_multiple_proteins<-matrix.mean.expression.ctrl.normalized[proteins,]
    mypalette <- colorRampPalette(c("white","black"))(n = 100)
    heatmap.2(matrix_heatmap_multiple_proteins,Rowv="FALSE",Colv="FALSE",dendrogram="none",scale="none",density.info="none",trace="none",key=TRUE, col= mypalette,cexRow = 0.8)
    
  }
  if(treatment == "rnase")
  {
    matrix_heatmap_multiple_proteins<-matrix.mean.expression.rnase.normalized[proteins,]
    mypalette <- colorRampPalette(c("white","black"))(n = 100)
    heatmap.2(matrix_heatmap_multiple_proteins,Rowv="FALSE",Colv="FALSE",dendrogram="none",scale="none",density.info="none",trace="none",key=TRUE, col= mypalette,cexRow = 0.8)
    
  }
  if(treatment == "both")
  {
    matrix_heatmap_multiple_proteins<-as.matrix(table.mean.expression.combined[table.mean.expression.combined$Gene == proteins[1],c(1:18)])
    for(i in 2:length(proteins))
    {
      matrix_heatmap_multiple_proteins<-rbind(matrix_heatmap_multiple_proteins,as.matrix(table.mean.expression.combined[table.mean.expression.combined$Gene == proteins[i],c(1:18)]))
    }
    mypalette <- colorRampPalette(c("white","black"))(n = 100)
    heatmap.2(matrix_heatmap_multiple_proteins,Rowv="FALSE",Colv="FALSE",dendrogram="none",scale="none",density.info="none",trace="none",key=TRUE, col= mypalette,cexRow = 0.7,cexCol=0.8)
    
  }
}





