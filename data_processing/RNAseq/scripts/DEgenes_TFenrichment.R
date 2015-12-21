

#Produces Figure 2B
#====

args=commandArgs(trailingOnly=TRUE)
DE_GENES=args[1]
DE_FILE=args[2]
TF_FILE=args[3]
FIG_OUT=args[4]
SIG=as.numeric(args[5])
print(args)

de_data=read.table(DE_FILE,header=TRUE)
dexp_genes=read.table(DE_GENES,header=TRUE)[,'gene']

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
pdf(FIG_OUT,height=6,width=3)
print(enrichmentPlot)
dev.off()








