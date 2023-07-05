## Re-analyze the Framingham Heart Study data considered by Demler et
## al. The FHS data is available by application from Natl Heart Lung
## and Blood Institute.

require(delong.alt)
source('misc.R')
fhs.cohort <- read.csv('./fhs/vr_riskscorepop_0.csv')
fhs.offspring <- read.csv('./fhs/vr_riskscorepop_1.csv')
colnames(fhs.cohort) <- colnames(fhs.offspring) <- tolower(colnames(fhs.cohort))
fhs <- rbind(fhs.cohort,fhs.offspring)
w <- subset(fhs,select=c(lnage,lntot,lnhdl,lnsbp,lndbp,diab126,smokes))
w <- as.matrix(w)
d <- fhs$cvd

alpha <- .05
glm0 <- glm(d~w-1,family=binomial('logit'))
summary(glm0)
significant.vars <- summary(glm0)$coefficients[,'Pr(>|z|)'] < alpha
w.full <- cbind(1,w)
glm0 <- glm(d~w.full-1,family=binomial('logit'))
summary(glm0)
significant.vars <- summary(glm0)$coefficients[,'Pr(>|z|)'] < alpha
CIs <- lapply(colnames(w.full)[significant.vars][-1], function(col) {
    w.red <- w[,colnames(w)!=col]
    h <- plogis
    h.1 <- function(x)h(x)*(1-h(x))
    h.2 <- function(x)h.1(x)*(1-2*h(x))
    out <- lapply(list(full=w.full,red=w.red), function(w) {
        x <- w[d==0,]; y <- w[d==1,]
        beta.hat <- coef(glm(d~w-1,family=binomial('logit')))
        obs <- auc.hat(x%*%beta.hat,y%*%beta.hat)
        params <- list(p=ncol(w),beta=beta.hat,link=h,link.deriv=h.1,link.deriv2=h.2)
        linearized <- lapply(list(delong= function(w,d,beta)rep(0,ncol(w)), proposed=NULL), function(deriv.fn)
            auc.index.linearize(w,d,beta.hat,infl.fn=function(w,d)infl.glm(w,d,params,terms.only=TRUE),deriv.fn=deriv.fn))
        linearized <- simplify2array(linearized)
        list(obs=obs,linearized=linearized)
    })
    diff.linearized <- out$full$linearized - out$red$linearized
    diff.obs <- out$full$obs - out$red$obs
    diff.var.hat <- apply(diff.linearized,2,function(iid)var(iid) / length(iid))
    c(diff.obs) / sqrt(diff.var.hat)
    CI.delong <- diff.obs + c(-1,1)*qnorm(1-alpha/2)*sqrt(diff.var.hat['delong'])
    CI.proposed <- diff.obs + c(-1,1)*qnorm(1-alpha/2)*sqrt(diff.var.hat['proposed'])
    rbind(delong=CI.delong,proposed=CI.proposed)
})
CIs <- abind::abind(CIs,along=3)
dimnames(CIs)[[3]] <- colnames(w.full)[significant.vars][-1]
CI.lengths <- apply(CIs,c(1,3),diff)


## produce figure for manuscript

extrafont::loadfonts()
png('fhs.png', width = 1024, height = 768, pointsize=15, family='CM Roman')
plot(0,type='n',xlim=range(CIs),ylim=c(0,7),xlab='difference in AUC',ylab='',yaxt='n')
axis(2, las=1, at=1:dim(CIs)[3], labels=colnames(w.full)[significant.vars][-1])
for(i in 1:dim(CIs)[3]) {
    m <- CIs[,,i]
    x0 <- m[,1]
    y0 <- i+c(0,1)/4
    x1 <- m[,2]
    segments(x0=x0,y0=y0,x1=x1,lty=1:2)
    tick.len <- .05
    segments(x0=x0,y0=y0-tick.len,y1=y0+tick.len)
    segments(x0=x1,y0=y0-tick.len,y1=y0+tick.len)
}
legend('topright',lty=1:2,legend=c('Delong','proposed'))
dev.off()



