mlmfit <-
function(jm1,data,Maxiter=200) {
  data_surv=data
  res1<-design_frailtyHL(jm1[[3]],data=data)
  kkk<-length(res1[[4]])
  if (length(res1[[4]])==1) res<-mlmfit1(jm1=jm1,data_surv,Maxiter=Maxiter)
  if (length(res1[[4]])==2) res<-mlmfit2(jm1=jm1,data_surv,Maxiter=Maxiter)
  return(res)
}
