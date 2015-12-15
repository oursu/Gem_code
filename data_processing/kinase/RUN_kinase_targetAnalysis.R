
#RMYPATH=/nfs/vendata/oursu/oana/GemPaper_2015-12-07/src/R/R-3.2.2/bin/R
#export R_LIBS=$R_LIBS:/nfs/vendata/oursu/oana/GemPaper_2015-12-07/src/R/my_libraries/
#installing packages: install.packages('knitr',lib='/nfs/vendata/oursu/oana/GemPaper_2015-12-07/src/R/my_libraries/',repos='http://cran.rstudio.com/')
#getting them: require(knitr)

args=commandArgs(trailingOnly=TRUE)
DATADIR=args[1]
require(knitr)
inrmd=paste(DATADIR,'/src/kinase/kinase_targetAnalysis.Rmd',sep='')
res=paste(DATADIR,'/results/data_processing/kinase/kinase_targetAnalysis',sep='')
knit2html(inrmd,res)
