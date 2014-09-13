d <- data.frame(lapply(read.csv('test.csv', sep=',', dec='.', stringsAsFactors=F),
                       as.numeric))

library(leaps)

x1 <- as.matrix(d[, 1:6])
y <- d$y

x2 <- as.matrix(d[, c(1, 1:6)])

x3 <- as.matrix(d[, c(1:6, 1)])


x <- x1
wt=rep(1,length(y))
force.in=NULL


force.out=NULL
intercept=FALSE
nvmax=8
nbest=1
warn.dep=TRUE
