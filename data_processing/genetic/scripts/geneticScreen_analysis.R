#Analysis of genetic hits
#========================================================

args=commandArgs(trailingOnly = TRUE)
screen_location=args[1]
pval_threshold=as.double(args[2])
SI_dev_threshold=as.double(args[3])
veh_survive=as.double(args[4])
geneticHitsFilename=args[5]
print(args)

#Read in data
data=read.table(screen_location,header=TRUE)

#GENETIC HITS
#Extract genetic hits: adj.P.Val, abs(SI=gem/veh-1)>=SI, veh survival
out_data=data
#remove controls
out_data=out_data[which(as.character(out_data$controlStatus)=='sample'),]

#SI=out_data$ratio
SI=(out_data$Dnorm.1+out_data$Dnorm.2)/(out_data$Vnorm.1+out_data$Vnorm.2)
SI_dev=SI-1
data.4=cbind(out_data,SI_dev=SI_dev)
survivors=which((data.4$Vnorm.1+data.4$Vnorm.2)/2>=veh_survive)

geneticHits=intersect(survivors,intersect(which(data.4$FDR<=pval_threshold),which(abs(data.4$SI_dev)>=SI_dev_threshold)))
genetic_hit_data=data.4[geneticHits,]
#Write down the genetic hits.
write.table(genetic_hit_data,file=paste(geneticHitsFilename,'.tab',sep=''),sep='\t',quote=F,row.names=F,col.names=T)
#Write down full analysis
write.table(data.4,file=paste(geneticHitsFilename,'.fullAnalysis_allGenes.tab',sep=''),sep='\t',quote=F,row.names=F,col.names=T)
print('done')
#=================
#END



