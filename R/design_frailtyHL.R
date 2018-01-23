design_frailtyHL <-
function (formula, data, weights, subset, na.action) 
{
    Call <- match.call()
    mc <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
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
    leftT<-0
    if (ncol(Y)==3) leftT<-1
    if (leftT==0) y <- matrix(Y[, 1], length(fr$Y), 1)
    if (leftT==0) L0 <- matrix(0, length(fr$Y), 1)
    if (leftT==1) y <- matrix(Y[, 2], length(fr$Y), 1)
    if (leftT==1) L0 <- matrix(Y[, 1], length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n <- nrow(x)
    p <- ncol(x)
    x1 <- x[1:n, 2:p]
    x2 <- matrix(x1, n, p - 1)
    x <- x2
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <- matrix(0, n, 1)
    if (leftT==0) del[, 1] <- censor <- Y[, 2]
    if (leftT==1) del[, 1] <- censor <- Y[, 3]
    SS <- FL$Subject
    res1 <- FrailtyMakeData(y, x, del, z, L0)
    y <- res1[1][[1]]
    x <- res1[2][[1]]
    del <- res1[3][[1]]
    z <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- matrix(0, p, 1)
    qcum <- cumsum(c(0, q))
    v_h <- matrix(0, qcum[nrand + 1], 1)
    res <- list(y, x, del, z, Mi, idx2, t2, di,beta_h,qcum,v_h,p,q,nrand,namesX)
    return(res)
}
