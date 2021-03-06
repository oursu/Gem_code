

#Produces Figure 2A, B, C
#====

args=commandArgs(trailingOnly=TRUE)
DATADIR=args[1]
FIG_OUT=args[2]
DE_GENES=args[3]
GENETIC_HITS=args[4]
EXPRESSED_KINASE_HITS=args[5]
DE_FILE=args[6]
FPKMmin=args[7]
TF_FILE=args[8]

print(args)
de_data=read.table(DE_FILE,header=TRUE)
print('read DE')
#Comparison with the genetic hits and kinase hits
#------

expressed_kinase_hits=read.table(EXPRESSED_KINASE_HITS,sep='\t')[,1]
#also read a table of the total genetic hits
genetic_hits=unique(setdiff(as.character(read.table(GENETIC_HITS,header=TRUE)[,1]),'NA'))
#also read a table of expressed genes (for filtering)
de_data=read.table(DE_FILE,header=TRUE)
expressed_genes=de_data[union(which(de_data$value_1>=FPKMmin),which(de_data$value_2>=FPKMmin)),'gene']
expressed_genetic_hits=intersect(genetic_hits,expressed_genes)
dexp_genes=read.table(DE_GENES,header=TRUE)[,'gene']

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

#plot for the figure
pdf(paste(FIG_OUT,'-Venn_Datasets.pdf',sep=''),width=4,height=4)
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

print('Now printing')                
print('Genetic and DE')
print(unique(intersect(expressed_genetic_hits,dexp_genes)))
print('Kinase and DE')
print(overlap_set_with_kinase(dexp_genes,expressed_kinase_hits))

tf_data_symbol=read.table(TF_FILE,header=FALSE)
colnames(tf_data_symbol)=c('TF','gene')
#items for the Fisher test
tfs_symbol=as.character(unique(tf_data_symbol$TF))
tfs_geneticHits=intersect(tfs_symbol,genetic_hits)
print('Genetic hit TFs')
print(tfs_geneticHits)


