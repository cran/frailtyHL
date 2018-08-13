frailty2.vs <-
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
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    
    del <- matrix(0, n, 1)
    del[, 1] <- censor <- Y[, 2]
    SS <- FL$Subject
    res1 <- FrailtyMakeData(y, x, del, z)
    y <- as.numeric(res1[1][[1]])
    x <- as.matrix(res1[2][[1]])
    del <- as.numeric(res1[3][[1]])
    z <- res1[4][[1]]

    	
    Mi <- as.matrix(res1[5][[1]])
    idx2 <- res1[6][[1]]
    t2 <- as.numeric(res1[7][[1]])
    di <- as.numeric(res1[8][[1]])



if(is.null(v)){
for (i in 1:nrand) v[[i]]<-rep(0, q[i])
}

if(is.null(alpha)){
alpha<-rep(0.1, nrand)
#alpha<-c(0.026, 1.036)
}


if(is.null(B)){
B<-as.numeric(frailtyHL(formula,data=data)$FixCoef[,1])
}


#print(z)

if(is.null(tun1)==FALSE|is.null(tun2)==FALSE){

if (model=="lognorm"){
result<-lognorm2.GCV(X=x,B,Z=z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del)
}

#if (model=="gamma"){
#result<-gammapv.GCV(X=x,B,Z=z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del)
#}
}

else{

if (model=="lognorm"){
result<-lognorm2.iter(X=x,B,Z=z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del)
}

#if (model=="gamma"){
#result<-gammapv.iter(X=x,B,Z=z,v,penalty,alpha,tun1,tun2,di,Mi,idx2,del)
#}
}


print("Result of variable selection in frailty model")
print("==Fitted model==")
print(paste("model :", model))
print(paste("penalty :", penalty))
print("formula :")
print(formula)
cn<-c("did not converge")
if(result$convergence==1){
cn<-c("converge")
}
print(cn)
print("==Fixed coefficients==")
print_mat<-cbind(result$beta,result$se)
rownames(print_mat)<-namesX
colnames(print_mat)<-c("Estimate","Std. Error")
print(print_mat)
print("==Dispersion parameter==")
print(result$alpha)
print("==Tuning parameter==")
print(result$tuning_parameter)
print("==BIC==")
print(as.numeric(result$BIC))


result$res1<-res1
return(result)

}
