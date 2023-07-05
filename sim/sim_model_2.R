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
    auc.f <- auc.full; auc.r <- auc.red # clean up
    p.r <- p.red; p.f <- p.full
    pi.1 <- 1-pi.0
    b <- 1 # clean up
    delta.0 <- .2 # clean up
    s.0.bar.r <- .1
    beta.star.r <- rep(b,p.r)
    beta.star.f <- rep(b,p.f)
    alpha <- qnorm(auc.f)/qnorm(auc.r)*sqrt(p.r/p.f)
    alpha.r <- 1/b/sqrt(p.r)*qnorm(auc.r)
    alpha.f <- 1/b/sqrt(p.f)*qnorm(auc.f)
    s.1.bar.r <- with(list(a=pi.1^2, b=2*pi.0*pi.1*s.0.bar.r-alpha.r^2, c=s.0.bar.r^2*pi.0^2-s.0.bar.r*alpha.r^2), 1/2/a*(-b + c(-1,1)*sqrt(b^2-4*a*c)))
    stopifnot(sum(s.1.bar.r>0)>0)
    s.1.bar.r <- s.1.bar.r[s.1.bar.r>0]
    gamma <- pi.0*p.r/p.f*s.0.bar.r+pi.0*delta.0+pi.1*p.r/p.f*s.1.bar.r
    delta.1 <- with(list(a=pi.1^2, b=2*pi.1*gamma-alpha.f^2, c=gamma^2-alpha.f^2*p.r/p.f*(s.0.bar.r+s.1.bar.r)-delta.0*alpha.f^2), 1/2/a*(-b + c(-1,1)*sqrt(b^2-4*a*c)))
    stopifnot(sum(delta.1>0)>0)
    delta.1 <- delta.1[delta.1>0]
    s.0.bar.f <- p.r/p.f*s.0.bar.r + delta.0
    s.1.bar.f <- p.r/p.f*s.1.bar.r + delta.1
    c(s.0.bar.r,s.1.bar.r)
    Sigma.0.r <- diag(epsilon+c(rep(0,p.r/2),rep(2*(s.0.bar.r-epsilon),p.r/2)))
    Sigma.1.r <- diag(epsilon+c(rep(0,p.r/2),rep(2*(s.1.bar.r-epsilon),p.r/2)))
    Sigma.0.f <- diag(c(diag(Sigma.0.r), rep(p.f*delta.0,p.f-p.r)))
    Sigma.1.f <- diag(c(diag(Sigma.1.r), rep(p.f*delta.1,p.f-p.r)))     
    sampler.red <- sampler.init.lda.2(p=p.r,Sigma.0=Sigma.0.r,Sigma.1=Sigma.1.r,beta.star=beta.star.r,pi.0=pi.0)
    sampler.full <- sampler.init.lda.2(p=p.f,Sigma.0=Sigma.0.f,Sigma.1=Sigma.1.f,beta.star=beta.star.f,pi.0=pi.0)
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
epsilon <- .2
auc.full <- .8
auc.red <- .77
terms <-  tryCatch(
    sim(n,p.full=p.full,auc.full=auc.full,auc.red=auc.red,epsilon=epsilon,pi.0=pi.0,reps=reps) ,
    error=function(e)NA)


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(terms,n,p.full,p.red,reps,auc.full,auc.red,pi.0,file=filename)



## Example shell command:
## for n in $(seq 1e3 1e3 8e3); do for reps in 300; do for pi0 in .05 .03; do /usr/bin/Rscript sim.R $n $reps 8 7 $pi0; done; done; done&


