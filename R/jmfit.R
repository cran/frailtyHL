jmfit <-
function(jm1,data1,jm2=NULL,data2=NULL,Maxiter) {
   if (is.null(jm2)==TRUE) {
    kkk<-length(jm1)
    if (kkk==2) res<-jmfit2(jm1[[1]],jm1[[2]],data1[[1]],data1[[2]],Maxiter)
    if (kkk==3) res<-jmfit3(jm1[[1]],jm1[[2]],jm1[[3]],data1[[1]],data1[[2]],Maxiter)
    } else {
      res<-jmfit_FM(jm1,data1,jm2,data2,Maxiter)
    }
    return(res)
}
