filelist <- dir()
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]


alpha <- .05
## file <- filelist[1] #!!!!!
coverage <- lapply(filelist, function(file) {
    ## print(file)
    load(file)
    out <- terms;     rm(terms) # clean up in sim.R
    if(identical(out,NA)) return(NA)
    lapply(split(out, out$stat), function(out.i) {
        stat <- unique(out.i$stat)
        obs <- out.i$obs
        true <- unname(unique(out.i$true))
        stopifnot(length(true)==1)
        var.hats <- as.matrix(out.i[,grep('^var',colnames(out.i))])
        test.stats <- (obs - true) / sqrt(var.hats)
        covered.sums <- colSums(abs(test.stats)<qnorm(1-alpha/2))
        length.sums <- colSums(2*sqrt(var.hats)*qnorm(1-alpha/2))
        data.frame( stat=stat, p.full=p.full,p.red=p.red,true=unique(true),pi.0=pi.0,n=n,estimator=sub('var.','',names(covered.sums)),reps=reps, covered.sums=covered.sums,length.sums=length.sums)
    })
    ## do.call(rbind, coverage)
})
mean(sapply(coverage,anyNA))
coverage <- coverage[!sapply(coverage,anyNA)]
coverage <- do.call(c, coverage)
## coverage.i <- split(coverage, names(coverage))[[1]]
coverage <- lapply(split(coverage, names(coverage)), function(coverage.i) {
    coverage.i <- do.call(rbind,coverage.i)
    coverage.i <- aggregate( cbind(reps,covered.sums,length.sums) ~ ., sum, data=coverage.i)
    rownames(coverage.i) <- NULL
    within(coverage.i, {
        fpr <- covered.sums/reps
        mean.length <- length.sums/reps
    })
})
## save.image('230619.RData')

out <- subset(coverage$diff,select=c(estimator,n,pi.0,fpr))
by.pi.0 <- lapply( split(out,out$pi.0), function(df) {
    out <- reshape(subset(df,select=-pi.0),direction='wide',timevar='estimator',idvar='n')
    colnames(out) <- sub('out.','',colnames(out))
    out
})
by.pi.0 <- rev(by.pi.0)
## op <- par(mfrow=rep(length(by.pi.0)/2,2))
require(extrafont)
loadfonts()
## png(paste0('coverage_diff.png'), width = 1024, height = 768, pointsize=15, family='CM Roman')
op <- par(mfrow=c(1,3), oma=c(2, 2, 0, 0)+0.1)
## out <- by.pi.0[[3]]
for(i in 1:length(by.pi.0)) {
    out <- by.pi.0[[i]]
    pi.1 <- 1-as.numeric(names(by.pi.0)[i])
    matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:4,col=1,xlab='',ylab='', main=bquote(pi[1] == .(pi.1)),ylim=c(.8,1))
    abline(h=1-alpha)
}
## expression(paste(pi[1],1-as.numeric(names(by.pi.0)[i])))
## main <- switch(estimand,red='AUC of reduced model', full='AUC of full model', diff='difference in AUCs')
## mtext(main, side=3, outer=TRUE)
mtext("sample size", side=1, outer=TRUE)
mtext("coverage", side=2, outer=TRUE)
estimator.names <- substr(colnames(out),5,nchar(colnames(out)))[-1]
estimator.names[estimator.names=='delong'] <- 'Delong'
legend('bottomright',legend=estimator.names,lty=1:4)
par(op)
## dev.off()


## all 3 estimands
for(estimand in c('red','full','diff')) {
    out <- subset(coverage[[estimand]],select=c(estimator,n,pi.0,fpr))
    by.pi.0 <- lapply( split(out,out$pi.0), function(df) {
        out <- reshape(subset(df,select=-pi.0),direction='wide',timevar='estimator',idvar='n')
        colnames(out) <- sub('out.','',colnames(out))
        out
    })
    by.pi.0 <- rev(by.pi.0)
    ## op <- par(mfrow=rep(length(by.pi.0)/2,2))
    png(paste0('coverage_',estimand,'.png'), width = 1024, height = 768, pointsize=15)
    op <- par(mfrow=c(1,3), oma=c(2, 2, 2, 0)+0.1)
    ## out <- by.pi.0[[3]]
    for(i in 1:length(by.pi.0)) {
        out <- by.pi.0[[i]]
        matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:4,col=1,xlab='',ylab='', main=1-as.numeric(names(by.pi.0)[i]),ylim=c(.8,1))
        abline(h=1-alpha)
    }
    main <- switch(estimand,red='AUC of reduced model', full='AUC of full model', diff='difference in AUCs')
    mtext(main, side=3, outer=TRUE)
    mtext("sample size", side=1, outer=TRUE)
    mtext("coverage", side=2, outer=TRUE)
    legend('bottomright',legend=substr(colnames(out),5,nchar(colnames(out)))[-1],lty=1:4)
    par(op)
    dev.off()
}



dd

## z.stats <- lapply(filelist, function(file) {
##     ## print(file)
##     load(file)
##     if(isTRUE(is.na(terms)))  return(NA)
##     terms <- as.data.frame(t(terms))
##     obs <- unlist(subset(terms,select=obs))
##     ## true <- unlist(subset(terms,select=true))
##     var.hats <- as.matrix(subset(terms,select=-c(obs,true)))
##     ## print(true)
##     ## test.stats <- (obs - true) / sqrt(var.hats)
##     ## n.covered <- colSums(abs(test.stats)<qnorm(1-alpha/2))
##     ## obs / sqrt(var.hats)
##     data.frame( p.full=p.full,p.red=p.red,adjust.size=round(adjust.size,3),n=n, (obs) / sqrt(var.hats) )
## })
## z.stats <- do.call(rbind,z.stats)
## ## print(sum(is.na(coverage)))
## std.error <- aggregate( cbind(delong,oracle,proposed) ~ ., sd, data=z.stats)
## std.error






## ## figure
## ## png('sim_logit.png')
## out <- subset(coverage,select=c(estimator,n,pi.0,fpr))
## by.adj.size <- lapply( split(out,out$adjust.size), function(df) {
##     out <- reshape(subset(df,select=-adjust.size),direction='wide',timevar='estimator',idvar='n')
##     colnames(out) <- sub('out.','',colnames(out))
##     out
## })
## out <- by.adj.size[[1]]
## matplot(out$n,subset(out,select=-n),pch=1,type='l',lty=1:3,col=1,xlab='n',ylab='coverage')
## legend('bottomright',legend=substr(colnames(out),5,nchar(colnames(out)))[-1],lty=1:3)
## abline(h=1-alpha)
## ## dev.off()

## table
require(xtable)
ns <- c(1e3,5e3,8e3)
fpr <- subset(coverage, n%in%ns, select=c(estimator,n,pi.0,fpr))
fpr <- xtabs(fpr ~ ., fpr)
fpr <- ftable(aperm(fpr,c(3,1,2)))
fpr <- round(fpr,3)
lengths <- subset(coverage, n%in%ns, select=c(estimator,n,pi.0,mean.length))
lengths <- xtabs(mean.length ~ ., lengths)
lengths <- ftable(aperm(lengths,c(3,1,2)))
lengths <- round(lengths,3)
out <- fpr
for(i in 1:ncol(fpr)) out[,i] <- paste0(format(fpr[,i],nsmall=3),' (',format(lengths[,i],nsmall=3),')')
names(attr(out,'row.vars'))[names(attr(out,'row.vars'))=='pi.0'] <- 'imbalance'
## save.file <- 'table_lda.tex'
sink(save.file)
xtableFtable(out)
sink()
lines <- scan(save.file,what='',sep='\n')
start.idx <- grep('begin\\{tabular\\}',lines)
end.idx <- grep('end\\{tabular\\}',lines)
## writeLines(lines[start.idx:end.idx],save.file)





