args=commandArgs(trailingOnly=TRUE)
inrmd=args[1]
res=args[2]
require(knitr)
system(paste('cd ',res,sep=''))
knit2html(inrmd)
print('Done')
