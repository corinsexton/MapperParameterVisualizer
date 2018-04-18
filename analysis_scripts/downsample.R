library(readr)

args <- commandArgs(TRUE)
filename <- args[1]

x <- read_tsv(filename, col_names = c('chr','pos','depth'))
x <- x[-1]

n <- 10000;
y <- aggregate(x,list(rep(1:(nrow(x)%/%n+1),each=n,len=nrow(x))),mean)[-1]
z<-y[order(y$depth),]

write_tsv(z,paste(filename, '.10000window',sep = ''))
