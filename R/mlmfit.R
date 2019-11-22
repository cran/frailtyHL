mlmfit <-
function(jm1,data, weights, subset, na.action,Maxiter=200) {
    Call <- match.call()
    mc <- match.call()
    indx <- match(c("data", "weights", "subset", "na.action"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) 
        stop("A formula argument is required")
  res1<-design_frailtyHL(jm1[[3]],data=data)
  kkk<-length(res1[[4]])
  data_surv<-data
  if (length(res1[[4]])==1) res<-mlmfit1(jm1=jm1,data_surv,Maxiter=Maxiter)
  if (length(res1[[4]])==2) res<-mlmfit2(jm1=jm1,data_surv,Maxiter=Maxiter)
  return(res)
}
