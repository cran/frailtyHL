jmfit <-
function(jm,data,jm2=NULL,data2=NULL,Maxiter) {
   if (is.null(jm2)==TRUE) {
    kkk<-length(jm)
    if (kkk==2) res<-jmfit2(jm[[1]],jm[[2]],data[[1]],data[[2]],Maxiter)
    if (kkk==3) res<-jmfit3(jm[[1]],jm[[2]],jm[[3]],data[[1]],data[[2]],Maxiter)
    } else {
      res<-jmfit_FM(jm,data,jm2,data2,Maxiter)
    }
    return(res)
}
