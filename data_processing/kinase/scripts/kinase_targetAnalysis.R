args=commandArgs(trailingOnly=TRUE)
DATADIR=args[1]
KINASE_HITS=args[2]
DE_FILE=args[3]
KINASE_THRESHOLD=args[4]
FPKMmin=args[5]
OUT=args[6]
TARGET_NAME_MAPPING=args[7]
print(args)

#Produces Figure 1C, 1D
#======================

require(reshape2)
require(pheatmap)

system(paste('mkdir',OUT))
KI_NAMES=paste(DATADIR,'/data/kinase/Table.S1.2.KI_names.txt',sep='')
KINASE_TARGETS=paste(DATADIR,'/data/kinase/Table.S1.3.KinaseTargets.txt',sep='')
GENETIC_ANALYSIS=paste(DATADIR,'/data/genetic/Table.S2.1.GeneticHitsProcessedData.txt.translatedR_org.Hs.eg.db_2011-08-07_entrez_to_geneSymbols',sep='')
FIG_OUT=paste(OUT,'/',basename(OUT),sep='')
KINASE_PHEN=paste(FIG_OUT,'.phen',sep='')

target_name_map=read.table(TARGET_NAME_MAPPING,header=TRUE,sep='\t')
rownames(target_name_map)=as.character(target_name_map[,'Name'])

kinase.data=read.table(KINASE_TARGETS,header=TRUE,sep='\t')
kinase.data=data.frame(kinase.data,target_kinase=target_name_map[as.character(kinase.data[,1]),'symbol'])

kinase.data.old=kinase.data
kinase.data=kinase.data[1,]
kinase.data=kinase.data[-1,]
nonvalue_cols=which(colnames(kinase.data.old) %in% c('compound.name.','target_kinase','siRNAscreen.SI'))

kinase.data=kinase.data[,-c(nonvalue_cols,ncol(kinase.data.old))]
for (item in unique(kinase.data.old$target_kinase)){
    rows=which(as.character(kinase.data.old$target_kinase)==as.character(item))
    keepvals=colMeans(kinase.data.old[rows,-nonvalue_cols])
    for (i in c(1:length(keepvals))){
    	keepvals[i]=max(0,keepvals[i])
    } 
    keep=data.frame(t(keepvals),target_kinase=item)
    kinase.data=rbind(kinase.data,keep)
}
rownames(kinase.data)=kinase.data$target_kinase

#Genetic hits files
genetic.data=read.table(GENETIC_ANALYSIS,header=TRUE)[,c('TraslatedGeneSymbol',
                                                         'ratio','FDR','controlStatus')]
#also read a table of expressed genes (for filtering)
de_data=read.table(DE_FILE,header=TRUE)
expressed_genes=de_data[union(which(de_data$value_1>=FPKMmin),which(de_data$value_2>=FPKMmin)),'gene']

#Convert kinase data to genetic names as well.
#Loop through kinase data and get corresponding genetic hit
kinase.data=data.frame(kinase.data,siRNAscreen.SI=NA)
unique_kinases=unique(kinase.data$target_kinase)
genetic.genes=as.character(genetic.data$TraslatedGeneSymbol)
for (kinase in unique_kinases){
  kinase_items=strsplit(kinase,'/')[[1]]
  correspondingGene=c()
  for (kinase_item in kinase_items){
    new_item=which(genetic.genes==toupper(kinase_item))
    if (length(new_item)>0){
      correspondingGene=c(correspondingGene,new_item) 
    }
  }
  if (length(correspondingGene)>0){
    kinase.data[which(as.character(kinase.data$target_kinase)==kinase),
                        'siRNAscreen.SI']=max(abs(1-as.numeric(as.character(genetic.data[correspondingGene,'ratio']))),na.rm=TRUE)
  }
}
#Targets of gemcitabine synergizers.
#----

#Get the synergizer drugs.
#========================
KI_dict=read.table(KI_NAMES,header=TRUE,sep='\t')
KI_dict[,'Ki_name']=gsub('/','.',gsub(' ','.',gsub(',','.',gsub('-','.',KI_dict[,'Ki_name']))))
rownames(KI_dict)=KI_dict[,'KI']
synergizer_ids=as.character(unique(read.table(KINASE_HITS,header=TRUE,sep='\t')[,'KI']))
synergizers=as.character(KI_dict[synergizer_ids,'Ki_name'])
synergy_cols=which(colnames(kinase.data) %in% synergizers)

#Make a heatmap of targets x synergizers. Note that for some pairs, data are missing. The data are missing from the original paper, it is not a bug in this code (http://www.nature.com/nbt/journal/v29/n11/full/nbt.2017.html)

kinh.target=kinase.data[,synergy_cols]

get_optimal_ordering=function(m){
  require(cba)
  if (dim(m)[1]==2){
   m.optimalRows=m
  }
  if (dim(m)[1]>2){
   d <- dist(as.matrix(m))
   hc <- hclust(d)
   co <- order.optimal(d, hc$merge)
   m.optimalRows=as.matrix(m)[co$order,]
  }
  if (dim(m)[2]==2){
   m.optimal=m.optimalRows
  }
  if (dim(m)[2]>2){
   d <- dist(as.matrix(t(m)))
   hc <- hclust(d)
   co <- order.optimal(d, hc$merge)
   m.optimal=m.optimalRows[,co$order]
  }
  return(m.optimal)
}
pheatmap(get_optimal_ordering(t(kinh.target)),cluster_cols=FALSE,cluster_rows=FALSE,
         show_colnames=FALSE)

#plot for the figure
pdf(paste(FIG_OUT,'KinaseHitHeatmap.pdf',sep=''),height=4,width=6)
pheatmap(get_optimal_ordering(t(kinh.target)),cluster_cols=FALSE,cluster_rows=FALSE,
         show_colnames=FALSE)
dev.off()

#Prepare a list of targets with associated values for the network input

for (DIRECTION in c('Down','UpAndDown')){

#first, count the targets per kinase inhibitor
if (DIRECTION=='Down'){
binarized.kinh.target=(kinh.target<=as.numeric(KINASE_THRESHOLD))
}
if (DIRECTION=='UpAndDown'){
binarized.kinh.target=(kinh.target<=as.numeric(KINASE_THRESHOLD))+(kinh.target>=(100+as.numeric(KINASE_THRESHOLD)))
}

write.table(binarized.kinh.target,file=paste(KINASE_PHEN,'-KinaseAffected',DIRECTION,'-BinarizedtargetsFile',sep=''),sep='\t',quote=F,row.names=T,col.names=T)
write.table(colSums(binarized.kinh.target,na.rm=TRUE),file=paste(KINASE_PHEN,'-KinaseAffected',DIRECTION,'-TotalTargetsPerKI',sep=''),sep='\t',quote=F,row.names=T,col.names=T)


phendata=data.frame(target=character(),value=double(),SI=double())
for (target in rownames(kinh.target)){
  si=kinase.data[target,'siRNAscreen.SI']
  if (DIRECTION=='Down'){
   val=min(kinh.target[target,],na.rm=TRUE)
   if (val<=as.numeric(KINASE_THRESHOLD)){
    phendata=rbind(phendata,data.frame(target=target,value=abs(100-val)/100,SI=si))
   }
  }
  if (DIRECTION=='UpAndDown'){
   val_min=min(kinh.target[target,],na.rm=TRUE)
   val_max=max(kinh.target[target,],na.rm=TRUE)
   val=max(abs(100-val_min),abs(100-val_max))
   if ((val>=as.numeric(KINASE_THRESHOLD))){
    phendata=rbind(phendata,data.frame(target=target,value=val/100,SI=si))
    }
  }
} 
write.table(phendata,file=paste(KINASE_PHEN,'-KinaseAffected',DIRECTION,sep=''),sep='\t',quote=F,row.names=F,col.names=F)

#Compare the kinase data with the genetic data
#===

expressed_kinase_hits=c()
unique_kinaseHits=unique(phendata[,'target'])
for (kinase in unique_kinaseHits){
  kinase_items=strsplit(kinase,'/')[[1]]
  for (kinase_item in kinase_items){
    if (kinase_item %in% expressed_genes){
      expressed_kinase_hits=c(expressed_kinase_hits,kinase)
    }
  }
}

print('expressed hits')
print(length(unique(expressed_kinase_hits)))
write.table(phendata[which(as.character(phendata$target) %in% expressed_kinase_hits),],
file=paste(KINASE_PHEN,'-KinaseAffected',DIRECTION,'-expressed',sep=''),sep='\t',quote=F,row.names=F,col.names=F)

#Keep only synergizer kinase inhibitors, get per gene best % activity reduced. Plot this best activity reduction vs abs(SI-1).
synergizer.aggregated=data.frame(siRNAscreen.SI=kinase.data[rownames(kinh.target),'siRNAscreen.SI'],target_kinase=kinase.data$target_kinase)

synergizer.aggregated=data.frame(synergizer.aggregated,max=0,min=0,percent_activity=0)
for (i in c(1:dim(synergizer.aggregated)[1])){
 target=synergizer.aggregated[i,'target_kinase']
 synergizer.aggregated[i,'max']=max(kinh.target[target,],na.rm=TRUE)
 synergizer.aggregated[i,'min']=min(kinh.target[target,],na.rm=TRUE)
 synergizer.aggregated[i,'percent_activity']=max(abs(100-synergizer.aggregated[i,'max']),abs(100-synergizer.aggregated[i,'min']))/100
}

print(synergizer.aggregated)
pdf(paste(FIG_OUT,'-KinaseAffected',DIRECTION,'KinaseVsGenetics.pdf',sep=''),height=6,width=6)
plot((100*synergizer.aggregated$percent_activity),
     synergizer.aggregated$siRNAscreen.SI,
     xlab='Percent activity change by kinase inhibitor',
  ylab='Genetic screen effect size',cex.lab=1.5,cex.axis=1.5)
abline(v=50,col='purple',lty=2)
abline(h=0.3,col='blue',lty=2)

plot((100*synergizer.aggregated$percent_activity),
          synergizer.aggregated$siRNAscreen.SI,
          xlab='Percent activity change by kinase inhibitor',
       ylab='Genetic screen effect size',cex.lab=1.5,cex.axis=1.5)
abline(v=50,col='purple',lty=2)
abline(h=0.3,col='blue',lty=2)
text((100*synergizer.aggregated$percent_activity),
               synergizer.aggregated$siRNAscreen.SI,label=synergizer.aggregated$target_kinase)
dev.off()

}
