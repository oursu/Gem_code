


#====

args=commandArgs(trailingOnly=TRUE)
DE_FILE=args[1]
SIG=as.numeric(args[2])
FPKMmin=as.numeric(args[3])
out=args[4]

#Produce a set of upregulated genes, and a set of downregulated genes 

#Get data
de_data=read.table(DE_FILE,header=TRUE)

#Function for picking DE genes
pick_DE_genes=function(qvalueThresh,foldChange,FPKMmin,FPKMdiff,data){
  
  OKFoldChange=which(abs(data$log2.fold_change.)>=log(foldChange,2))
  OKFPKMmin=union(which(data$value_1>=FPKMmin),which(data$value_2>=FPKMmin))
  OKFPKMdiff=which(abs(data$value_1-data$value_2)>=FPKMdiff)
  OKqvalue=which(data$q_value<=qvalueThresh)
  
  #For DE, keep only OK genes
  OKtestgenes=which(as.character(data$status)=='OK')
  DEidx=intersect(OKqvalue,intersect(OKFPKMdiff,intersect(OKFPKMmin,OKFoldChange)))
  
  #Hits:
  DEgene_data=data[DEidx,]

  #split into up and downreg
  columns_of_interest=c('gene','log2.fold_change.',
                        'value_1','value_2')
  deGenes=list()
  deGenes[['up']]=DEgene_data[which(as.numeric(as.character(DEgene_data[,'log2.fold_change.']))>0),columns_of_interest]
  deGenes[['down']]=DEgene_data[which(as.numeric(as.character(DEgene_data[,'log2.fold_change.']))<0),columns_of_interest]
  return(deGenes)
}

#Get DE genes
DE_genes=pick_DE_genes(SIG,1.5,FPKMmin,1,de_data)
print(DE_genes)
write.table(DE_genes[['up']],file=paste(out,'.UpReg.txt',sep=''),quote=F,row.names=F,col.names=T,sep='\t')
write.table(DE_genes[['down']],file=paste(out,'.DownReg.txt',sep=''),quote=F,row.names=F,col.names=T,sep='\t')
write.table(rbind(DE_genes[['down']],DE_genes[['up']]),file=paste(out,'.UpAndDownReg.txt',sep=''),quote=F,row.names=F,col.names=T,sep='\t')


