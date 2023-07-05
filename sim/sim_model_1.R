args <- commandArgs(trailingOnly=TRUE)
n <- as.numeric(args[1])
reps <- as.numeric(args[2])
p.full <- as.numeric(args[3])
p.red <- as.numeric(args[4])
pi.0 <- as.numeric(args[5])

suppressPackageStartupMessages({
    require(mvtnorm)
    require(parallel)
    require(numDeriv)
})
start <- Sys.time()
source('../../misc.R')
sim <- function(n,p.full,auc.full,auc.red,epsilon=NULL,pi.0=NULL,adjust.size=NULL,reps=3e2,bootstrap.reps=1e2) {
    p.red <- p.full - 1
    pi.1 <- 1-pi.0
    Sigma.0.red <- Sigma.1.red <- diag(p.red)
    Sigma.0.full <- diag(c(rep(1,p.red),epsilon))
    Sigma.1.full <- diag(p.full)
    b.red <- sqrt(2/p.red)*qnorm(auc.red)
    b.full2 <- with(list(a=(pi.0*epsilon+pi.1)^2,
                      b=2*b.red^2*p.red*(pi.0*epsilon+pi.1) - (epsilon+1)*qnorm(auc.full)^2,
                      c=b.red^4*p.red^2-qnorm(auc.full)^2*2*b.red^2*p.red),
                 1/(2*a) * (-b + c(-1,1)*sqrt(b^2-4*a*c)))
    b.full2 <- min(b.full2[!(b.full2<0)])
    b.full <- sqrt(b.full2)
    beta.red <- rep(b.red,p.red)
    beta.full <- c(beta.red,b.full)
    sampler.red <- sampler.init.lda.2(p=p.red,Sigma.0=Sigma.0.red,Sigma.1=Sigma.1.red,beta.star=beta.red,pi.0=pi.0)
    sampler.full <- sampler.init.lda.2(p=p.full,Sigma.0=Sigma.0.full,Sigma.1=Sigma.1.full,beta.star=beta.full,pi.0=pi.0)
    dataset.full <- sampler.full$params; dataset.red <- sampler.red$params
    terms <- replicate(reps, simplify=FALSE, expr={
        xd <- sampler.full$sample(n)
        dataset.full$x <- xd$x; dataset.full$d <- xd$d 
        dataset.red$x <- xd$x[,1:p.red]; dataset.red$d <- xd$d
        paired <- lapply(list(full=dataset.full,red=dataset.red), function(dataset) {
            with(dataset, {
                x.0 <- x[d==0,]; x.1 <- x[d==1,]
                beta.hat <- coefs.lda(x.0,x.1)
                obs <- auc.hat(x.0%*%beta.hat,x.1%*%beta.hat) 
                linearized <- lapply(list(delong= function(x,d,beta)rep(0,ncol(x)),oracle=function(x,d,beta)deriv.star, proposed=NULL), function(deriv.fn)
                    auc.index.linearize(x,d,beta.hat,infl.fn=function(x,d)infl.lda(x,d,params=NULL,var.equal=FALSE),deriv.fn=deriv.fn))
                linearized <- simplify2array(linearized)
                list(obs=obs,true=auc,linearized=linearized)
            })
        })
        auc.full <- paired$full
        auc.red <- paired$red
        auc.diff <- mapply(`-`,auc.full,auc.red)
        out <- sapply(list(full=auc.full,red=auc.red,diff=auc.diff), function(auc) {
            var.hat <- apply(auc$linearized,2,function(iid)var(iid) / length(iid))
            c(obs=auc$obs,true=auc$true,var=var.hat)
        })
        out <- as.data.frame(t(out))
        out$stat <- rownames(out)
        auc.bs <- replicate(bootstrap.reps, {           
            idx <- sample(1:nrow(dataset.full$x),replace=TRUE)
            paired <- lapply(list(full=dataset.full,red=dataset.red), function(dataset) {
                with(dataset, {
                    x.bs <- x[idx,]; d.bs <- d[idx]
                    if(mean(d.bs)*mean(1-d.bs)==0) {
                        warning('bootstrap sample all cases or all controls')
                        return(NA)
                        }
                    x.0.bs <- x[d.bs==0,,drop=FALSE]; x.1.bs <- x.bs[d.bs==1,,drop=FALSE]
                    beta.hat.bs <- coefs.lda(x.0.bs,x.1.bs)
                    auc.hat(x.0.bs%*%beta.hat.bs,x.1.bs%*%beta.hat.bs)
                })
            })
            c(full=paired$full,red=paired$red)
        })
        auc.bs <- rbind(auc.bs,diff=auc.bs['full',]-auc.bs['red',])
        var.bootstrap <- apply(auc.bs,1,var,na.rm=TRUE)
        var.bootstrap <- as.data.frame(var.bootstrap)
        var.bootstrap$stat <- rownames(var.bootstrap)
        out <- merge(out,var.bootstrap)
    })
    do.call(rbind,terms)
}
epsilon <- .00
auc.full <- .8
auc.red <- .77
terms <-  tryCatch(
    sim(n,p.full=p.full,auc.full=auc.full,auc.red=auc.red,epsilon=epsilon,pi.0=pi.0,reps=reps) ,
    error=function(e)NA)




filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(terms,n,p.full,p.red,reps,auc.full,auc.red,pi.0,file=filename)



## example shell call:
## for n in $(seq 1e3 1e3 8e3); do for reps in 300; do for pi0 in .05 .03; do /usr/bin/Rscript sim.R $n $reps 8 7 $pi0; done; done; done&


