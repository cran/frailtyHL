design_hglm <-
function(formula, DataMain, BinomialDen=NULL,weights, subset, na.action) {
    mc <- match.call()
    formulaMean<-formula
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   qcum <- cumsum(c(0, q))
   v_h<-matrix(0,qcum[nrand+1],1)
   res<-list(y, x,z,v_h,n,p,q,namesX)
}
