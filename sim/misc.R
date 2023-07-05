


## sample LDA data under first model discussed in simulation section
sampler.init.lda.1 <- function(p,s,d=NULL,auc,pi.0,adjustment.size=NULL) {
    pi.1 <- 1-pi.0
    if(is.null(d))
        d <- sqrt(1+mean(s))/sqrt(p)/(pi.1+pi.0*mean(s)) * qnorm(auc)
    beta.star <- rep(1,p)*d
    Sigma.1 <- diag(p)
    Sigma.0 <- diag(s)
    Sigma.pi <- pi.0*Sigma.0+pi.1*Sigma.1
    Sigma <- Sigma.0+Sigma.1
    mu.0 <- rep(0,p)
    mu.1 <- as.numeric((pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star)
    if(is.null(adjustment.size)) {
        adjustment.vec <- dnorm(d*sqrt(p)*(pi.1+pi.0*mean(s))/sqrt(1+mean(s))) * (pi.1-pi.0)/sqrt(p)/(1+mean(s))^(3/2)  *  sqrt(s/pi.0+1/pi.1)/(pi.1+pi.0*s)*(s-mean(s))
        adjustment.size <- sqrt(sum(adjustment.vec^2))
        }
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0,beta=beta.star,auc=auc,pi.0=pi.0,adjustment.size=adjustment.size)
    params$deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    sample <- function(n) {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        if(n.0*n.1==0)warning(call.=TRUE)
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        d <- c(rep(0,n.0),rep(1,n.1))
        return(list(x=x,d=d,n.0=n.0,n.1=n.1))
    }
    return(list(params=params,sample=sample))
}




## sample LDA data under second model discussed in simulation section
sampler.init.lda.2 <- function(p,Sigma.0,Sigma.1,beta.star,auc=NULL,pi.0,adjustment.size=NULL) {
    pi.1 <- 1-pi.0
    mu.0 <- rep(0,p)
    mu.1 <- as.numeric((pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star)
    if(is.null(auc)) auc <- pnorm( sum(beta.star^2*(pi.0*diag(Sigma.0)+pi.1*diag(Sigma.1))) / sqrt(sum(beta.star^2*(diag(Sigma.0+Sigma.1)))) )
    if(is.null(adjustment.size)) {
        quad.0 <- as.numeric(t(beta.star)%*%Sigma.0%*%beta.star)
        quad.1 <- as.numeric(t(beta.star)%*%Sigma.1%*%beta.star)
        quad <- as.numeric(t(beta.star)%*%(Sigma.0/2+Sigma.1/2)%*%beta.star)
        quad.pi <- as.numeric(t(beta.star)%*%(pi.0*Sigma.0+pi.1*Sigma.1)%*%beta.star)
        adjustment.vec <- dnorm(quad.pi/sqrt(2*quad)) * (pi.0-pi.1)/(2*quad)^(3/2) * (  sqrt(diag(Sigma.0)/pi.0+diag(Sigma.1)/pi.1) * (quad.0*diag(Sigma.1)-quad.1*diag(Sigma.0)) / (pi.0*diag(Sigma.0)+pi.1*diag(Sigma.1)) * beta.star  )
        adjustment.size <- sqrt(sum(adjustment.vec^2))
        }    
    params <- list(mu.0=mu.0,mu.1=mu.1,Sigma.0=Sigma.0,Sigma.1=Sigma.1,pi.0=pi.0,beta=beta.star,auc=auc,pi.0=pi.0,adjustment.size=adjustment.size)
    params$deriv.star <- auc.deriv.lda.gaussian(beta.star,params)
    sample <- function(n) {
        n.0 <- rbinom(1,n,pi.0); n.1 <- n-n.0
        if(n.0*n.1==0)warning(call.=TRUE)
        x.0 <- rmvnorm(n.0,mu.0,Sigma.0)
        x.1 <- rmvnorm(n.1,mu.1,Sigma.1)
        x <- rbind(x.0,x.1)
        d <- c(rep(0,n.0),rep(1,n.1))
        return(list(x=x,d=d,n.0=n.0,n.1=n.1))
    }
    return(list(params=params,sample=sample))
}


## influence function for heteroscedastic LDA model
infl.lda <- function(x,d,params,var.equal=TRUE,terms.only=TRUE) {
    g <- d # clean up
    n <- length(g) # clean up redundancies here
    n.1 <- sum(g); n.0 <- n-n.1
    mu.0 <- params$mu.0
    mu.1 <- params$mu.1
    pi.0 <- params$pi.0
    x.0 <- x[g==0,]; x.1 <- x[g==1,]
    if(is.null(mu.0))mu.0 <- colMeans(x.0)
    if(is.null(mu.1))mu.1 <- colMeans(x.1)
    if(is.null(pi.0))pi.0 <- n.0/n
    pi.1 <- 1-pi.0
    if(var.equal) {
        Sigma <- params$Sigma
        if(is.null(Sigma))Sigma <- (var(x.0)*n.0+var(x.1)*n.1)/n
    } else {
        Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
        if(is.null(Sigma.0))Sigma.0 <- var(x.0)
        if(is.null(Sigma.1))Sigma.1 <- var(x.1)
        Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
    }
    x <- t(x)
    g <- sort(g)
    infl <- t(g*t(solve(Sigma)%*%(x-mu.1)*n/n.1)  -  (1-g)*t(solve(Sigma)%*%(x-mu.0)*n/n.0))
    if(terms.only) return(infl) else return(rowMeans(infl))
}


## influence function for heteroscedastic GLM estimates
infl.glm <- function(x,g,params,terms.only=TRUE) {
    ## browser()
    p <- params$p
    score <- function(x,g,params) { # x in model matrix format
        h <- params$link; h.1 <- params$link.deriv
        eta <- as.numeric(x%*%params$beta)
        t( (g-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
    }
    fi <- function(x,g,params) { # Fisher inf, x in model matrix format
        h <- params$link; h.1 <- params$link.deriv; h.2 <- params$link.deriv2
        n <- nrow(x)
        eta <- as.numeric(x%*%params$beta)
        denom <- h(eta)*(1-h(eta))
        -t((h.2(eta)*(g-h(eta))/denom - h.1(eta)/denom^2 * (g*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x / n
    }
    terms <- solve(fi(x,g,params))%*%score(x,g,params)
    if (terms.only) return(terms) else return(rowMeans(terms))
}


## compute LDA coefficient estimates
coefs.lda <- function(x.0,x.1,params=NULL) {
    n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
    mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
    Sigma <- params$Sigma
    if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
    pi.0 <- params$pi.0 ## for heteroskedasticity
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(pi.0))pi.0 <- 1/2 
    if(is.null(Sigma.0))Sigma.0 <- Sigma
    if(is.null(Sigma.1))Sigma.1 <- Sigma
    Sigma <- pi.0*Sigma.0 + (1-pi.0)*Sigma.1
    Sigma.inv <- solve(Sigma)
    beta.hat <- Sigma.inv%*%(mu.1.hat-mu.0.hat)
    return(beta.hat)
}

## compute observed AUC
auc.hat <- function(x,y,ties=FALSE)
    if(!ties) {
        if(as.numeric(length(x))*length(y)<1e8)
            return(mean(outer(x,y,'<')))
        else {
                total <- 0
                for(i in 1:length(x))total <- total + sum(x[i]<y)
                return(total/length(x)/length(y))
            }
    } else return(NA) ## todo

