

#Produces Figure 2A
#====

args=commandArgs(trailingOnly=TRUE)
GOfile=args[1]
OUT=args[2]
SIG=args[3]
print(args)

#GO enrichment for differentially expressed genes
#--------

#Use DAVID website to perform GO enrichment analysis for Biological Process. Analysis performed on 2015-03-06. Below are the results.

require(ggplot2)
display_GO=function(go_data_name,out,h,w,SIG){
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

display_GO(GOfile,OUT,3,6,SIG)



