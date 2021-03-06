

#Produces Figure 2A, B, C
#====

args=commandArgs(trailingOnly=TRUE)
DATADIR=args[1]
FIG_OUT=args[2]
DE_FILE=args[3]
GENETIC_HITS=args[4]
EXPRESSED_KINASE_HITS=args[5]
SIG=as.numeric(args[6])
FPKMmin=as.numeric(args[7])
TF_FILE=args[8]

#2. Analysis
#-------
#Produce a set of upregulated genes, and a set of downregulated genes for GO enrichment analysis.

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
out=FIG_OUT
write.table(DE_genes[['up']],file=paste(out,'.UpReg.txt',sep=''),quote=F,row.names=F,col.names=T,sep='\t')
write.table(DE_genes[['down']],file=paste(out,'.DownReg.txt',sep=''),quote=F,row.names=F,col.names=T,sep='\t')

#GO enrichment for differentially expressed genes
#--------

#Use DAVID website to perform GO enrichment analysis for Biological Process. Analysis performed on 2015-03-06. Below are the results.

display_GO=function(go_data_name,out,h,w){
  #go_data_name=paste(out,'DownReg.GO.txt',sep='')
  go=read.table(go_data_name,header=TRUE,sep='\t')
  sig.go=go[which(as.numeric(as.character(go$Benjamini))<=SIG),]
  sig.go=data.frame(sig.go,logp=log(sig.go$Benjamini,10))
  sig.go=data.frame(sig.go,GOTerm=factor(sig.go$Term,levels=sig.go[order(sig.go$logp,decreasing=TRUE),'Term']))
  require(ggplot2)
  go_plot=ggplot(sig.go,aes(x=GOTerm,y=-logp))+geom_bar(stat='identity')+coord_flip()+theme_bw() 
  print(go_plot)
  print(sig.go)
  pdf(out,height=h,width=w)
  print(go_plot)
  dev.off()
}

#GO enrichment for downregulated genes

display_GO(paste(out,'.DownReg.GO.txt',sep=''),
           paste(FIG_OUT,'-DownRegGO.pdf',sep=''),
           3,6)

#GO enrichment for upregulated genes

display_GO(paste(out,'.UpReg.GO.txt',sep=''),
           paste(FIG_OUT,'-UpRegGO.pdf',sep=''),
           3,6)

print('up')
print(dim(DE_genes[['up']])[1])
print('down')
print(dim(DE_genes[['down']])[1])

#Comparison with the genetic hits and kinase hits
#------

expressed_kinase_hits=read.table(EXPRESSED_KINASE_HITS,sep='\t')[,1]
#also read a table of the total genetic hits
genetic_hits=unique(setdiff(as.character(read.table(GENETIC_HITS,header=TRUE)[,1]),'NA'))
#also read a table of expressed genes (for filtering)
de_data=read.table(DE_FILE,header=TRUE)
expressed_genes=de_data[union(which(de_data$value_1>=FPKMmin),which(de_data$value_2>=FPKMmin)),'gene']
expressed_genetic_hits=intersect(genetic_hits,expressed_genes)
dexp_genes=c(as.character(DE_genes[['up']][,'gene']),
           as.character(DE_genes[['down']][,'gene']))
overlap_set_with_kinase=function(setGenes,kinase_dataset){
  counted=c()
  for (kinase in kinase_dataset){
    kinase_items=strsplit(kinase,'/')[[1]]
    for (kinase_item in kinase_items){
      if (kinase_item %in% setGenes){
        counted=c(counted,kinase)
      }
    }
  }
  return(unique(counted))
}
#ready for the comparison
comparison=data.frame(genetic=c(0,0,0),kinase=c(0,0,0),DE=c(0,0,0))
rownames(comparison)=c('genetic','kinase','DE')
comparison['genetic','genetic']=length(expressed_genetic_hits)
comparison['genetic','kinase']=comparison['kinase','genetic']=length(overlap_set_with_kinase(expressed_genetic_hits,expressed_kinase_hits))
comparison['genetic','DE']=comparison['DE','genetic']=length(unique(intersect(expressed_genetic_hits,dexp_genes)))
comparison['kinase','kinase']=length(unique(expressed_kinase_hits))
comparison['kinase','DE']=comparison['DE','kinase']=length(overlap_set_with_kinase(dexp_genes,expressed_kinase_hits))
comparison['DE','DE']=length(unique(dexp_genes))
everything=length(overlap_set_with_kinase(expressed_genetic_hits,overlap_set_with_kinase(dexp_genes,expressed_kinase_hits)))

#And make a nice Venn
require(VennDiagram)
grid.newpage()
draw.triple.venn(comparison['genetic','genetic'],
                 comparison['kinase','kinase'],
                 comparison['DE','DE'], 
                 comparison['genetic','kinase'], 
                 comparison['kinase','DE'], 
                 comparison['genetic','DE'], 
                 everything,category=c('Genetic','Kinase','Differentially expressed'),
                 fill=c('lightblue','Purple','orange'))

#plot for the figure
pdf(paste(FIG_OUT,'-Venn_Datasets.pdf',sep=''))
grid.newpage()
draw.triple.venn(comparison['genetic','genetic'],
                 comparison['kinase','kinase'],
                 comparison['DE','DE'], 
                 comparison['genetic','kinase'], 
                 comparison['kinase','DE'], 
                 comparison['genetic','DE'], 
                 everything,category=c('Genetic','Kinase','Differentially expressed'),
                 fill=c('lightblue','Purple','orange'))
dev.off()
                
print('Genetic and DE')
print(unique(intersect(expressed_genetic_hits,dexp_genes)))
print('Kinase and DE')
print(overlap_set_with_kinase(dexp_genes,expressed_kinase_hits))

#TFs enriched as regulators of the differentially expressed genes
#--------------

#Universe: all genes in the differential expression analysis. 
#Contingency table: gene is diffExpr, gene is target of TF.


tf_data=read.table(TF_FILE,header=FALSE)
colnames(tf_data)=c('TF','gene')
#items for the Fisher test
tfs=as.character(unique(tf_data$TF))
all_genes=as.character(unique(de_data$gene))

tf_fisher=function(tf,all_genes,dexp_genes,tf_data){
    TFgenes=as.character(unique(tf_data[which(as.character(tf_data$TF)==tf),'gene']))
    m=data.frame(DE=c(0,0),notDE=c(0,0))
    rownames(m)=c('TF','notTF')
    #fill matrix
    m['TF','DE']=length(intersect(dexp_genes,TFgenes))
    m['notTF','DE']=length(setdiff(dexp_genes,TFgenes))
    m['TF','notDE']=length(setdiff(TFgenes,dexp_genes))
    m['notTF','notDE']=length(setdiff(all_genes,union(TFgenes,dexp_genes)))
    test=fisher.test(m)
    return(data.frame(TF=tf,
                      TF_DE=m['TF','DE'],
                      notTF_DE=m['notTF','DE'],
                      TF_notDE=m['TF','notDE'],
                      notTF_notDE=m['notTF','notDE'],
                      p=test[['p.value']],
                      OR=test[['estimate']],
                      conf.low=test[['conf.int']][1],
                      conf.high=test[['conf.int']][2]))   
}

TFenrichments=data.frame(TF='TF',TF_DE=0,notTF_DE=0,
                         TF_notDE=0,notTF_notDE=0,
                         p=0,OR=0,conf.low=0,conf.high=0)
TFenrichments=TFenrichments[-1,]
for (tf in tfs){
  TFenrichments=rbind(TFenrichments,tf_fisher(tf,all_genes,dexp_genes,tf_data))
}
TFenrichments=data.frame(TFenrichments,BH=p.adjust(TFenrichments[,'p'],method='BH'),
                         N=TFenrichments$TF_DE)
N_THRESHOLD=100
TFenrichments=TFenrichments[which(as.numeric(as.character(TFenrichments$N))>=N_THRESHOLD),]
TFenrichments=data.frame(TFenrichments,TF_withN=paste(TFenrichments$TF,' | N=',
                                                      TFenrichments$N))

#plot the enrichments, color by significance, rank by enrichment value
significant=rep('Not significant',times=dim(TFenrichments)[1])
significant[intersect(which(TFenrichments$BH<=SIG),
                      which(TFenrichments$conf.low>1))]='Significant'

TFenrichments=data.frame(TFenrichments,rankedTF=factor(TFenrichments$TF_withN,levels=TFenrichments[order(TFenrichments$OR),'TF_withN']),sig=significant)
toptfs=order(TFenrichments$OR,decreasing=TRUE)[1:50]
TFenrichments=TFenrichments[toptfs,]
require(ggplot2)
enrichmentPlot=ggplot(TFenrichments,aes(x=rankedTF,y=OR))+geom_bar(stat = "identity",fill="gray", colour="black") + coord_flip()+theme_bw()+ylim(0,3)+theme(axis.text.y  = element_text(size=6))+xlab('Transcription factor')+ylab('Odds ratio')
print(enrichmentPlot)
pdf(paste(FIG_OUT,'-TFenrichmentDEgenes.pdf',sep=''),height=6,width=3)
print(enrichmentPlot)
dev.off()


#Find the TFs that are genetic hits.

tf_data_symbol=tf_data
colnames(tf_data_symbol)=c('TF','gene')
#items for the Fisher test
tfs_symbol=as.character(unique(tf_data_symbol$TF))
tfs_geneticHits=intersect(tfs_symbol,genetic_hits)

print('genetic hit tfs')
print(tfs_geneticHits)



