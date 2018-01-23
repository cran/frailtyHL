frailty.vs <-
function(formula,model,penalty,data,B=NULL,v=NULL,alpha=NULL,tun1=NULL,tun2=NULL){

    Call <- match.call()
    mc <- match.call()
    indx <- match(c("formula", "data"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) 
        stop("A formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(subbars(formula), special)
    m <- eval(temp)
    Terms <- attr(m, "terms")
    Y <- model.extract(m, "response")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(formula, special)
    Terms <- temp[[2]]
    formula1 <- paste(paste(Terms[[2]][[2]], Terms[[3]], sep = "~")[[2]], 
        paste(Terms[[3]])[3], sep = "+")
    formula1 <- formula(formula1)
    fr <- FrailtyFrames(mc, formula1, contrasts)
    namesX <- names(fr$fixef)
    namesX <- namesX[-1]
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formula1, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(Y[, 1], length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n <- nrow(x)
    p <- ncol(x)
    x1 <- x[1:n, 2:p]
    x2 <- matrix(x1, n, p - 1)
    x <- x2
    n <- nrow(x)
    p <- ncol(x)
    alphas<-alpha
    nrand <- length(z)
    n.alpha<-length(alphas)
    if (n.alpha<=nrand) result<-frailty1.vs(formula=formula,model=model,penalty=penalty,data=data,B=B,v=v,alpha=alpha,tun1=tun1,tun2=tun2)
    else result<-frailty2.vs(formula=formula,model=model,penalty=penalty,data=data,B=B,v=v,alpha=alpha,tun1=tun1,tun2=tun2)
    return(result)
}
