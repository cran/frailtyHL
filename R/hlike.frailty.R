hlike.frailty <-
function(formula, data, inits, order=1, frailty.cov="none",
subHazard=FALSE, alpha=.05, MAX.ITER=100, TOL=1E-6) {
# Create a sorted data matrix from the given formula and dataset.
# -----------------------------------------------------------------------------
mf <- model.frame(formula, data)
if(!is(model.response(mf), "CmpRsk") | is.na(pmatch("cluster", names(mf))))
stop("Invalid formula")
resp <- pmatch("CmpRsk", names(mf))
clust <- pmatch("cluster", names(mf))
pred <- names(mf)[setdiff(1:ncol(mf), c(resp, clust))]
temp.data <- data.frame(unclass(mf[,resp]), mf[,pred], mf[,clust])
names(temp.data) <- c("time", "index", pred, "cluster")
sort.data <- temp.data[order(temp.data$time),]
N <- nrow(sort.data) # Total sample size
n <- length(unique(sort.data$cluster)) # Number of clusters
m <- max(sort.data$index) # Number of event types
p <- length(pred) # Number of predictors
# Create design matrices for fixed covariates and the frailty term
X <- as.matrix(sort.data[,pred]); colnames(X)=pred
Z <- model.matrix(~ 0 + factor(sort.data$cluster)); colnames(Z)=NULL
# Number of frailty parameters and derivative of the covariance matrix (dSigma)
# -----------------------------------------------------------------------------
if(m==1 & frailty.cov!="none") stop("Invalid frailty.cov argument")
if(subHazard & frailty.cov!="none") stop("Invalid frailty.cov argument")
s <- switch(frailty.cov,
none = 1,
independent = m,
exchangeable = m+1,
unstructured = sum(1:m),
stop("Invalid frailty.cov argument"))
dSigma <- function(q) {
dSigma <- matrix(0, nrow=min(m,s), ncol=min(m,s))
zero=rep(0,s); zero[q] <- 1
if(frailty.cov=="independent") { dSigma = diag(zero)
} else if(frailty.cov=="exchangeable") {
if(q<=nrow(dSigma)) dSigma <- diag(zero[-s])
if(q>nrow(dSigma)) dSigma[lower.tri(dSigma) | upper.tri(dSigma)]=1
} else if(frailty.cov=="unstructured") {
dSigma[lower.tri(dSigma, diag=TRUE)] <- zero
dSigma[upper.tri(dSigma)] <- dSigma[lower.tri(dSigma)]
} else if(frailty.cov=="none") dSigma = 1
return(dSigma)
}


# If the subhazard model is being fitted (subHazard==TRUE) then the event of
# interest is event 1. The function G is the Kaplan-Meier estimate of the
# survival function for the censoring distribution.
# -----------------------------------------------------------------------------
if(subHazard) {
m = 1
fit.censor <- survfit(Surv(sort.data$time, sort.data$index==0)~1)
G <- function(t) {
# The times argument of summary.survfit requires that t be sorted.
# After getting survival times at each t, reorder s to match the
# original order of t. Return a vector of 1s if there is no censoring.
s <- summary(fit.censor, times=sort(t))$surv
if(any(sort.data$index==0)) { return(s[match(t, sort(t))])
} else return(rep(1, length(t)))
}
}
# Unique event times (y), number of events per unique time (d), event indicator
# (delta) and risk matrix (atRisk) for event type k, k=1,...,m
#
# The atRisk matrix is a matrix of the risk set. Columns are unique event times
# with one row for each observation. A 1 denotes that the subject is still at
# risk; with sorted data resembles a lower-triangular matrix.
# -----------------------------------------------------------------------------
# Calculate, y, d, delta, and atRisk for each event type
prepData <- function(event) {
y <- unique(with(sort.data, time[index==event]))
d <- with(sort.data,tapply(index[index==event],time[index==event],length))
num.y <- length(y)
delta <- ifelse(sort.data$index==event,1,0)
atRisk <- matrix(0, nrow=N, ncol=num.y)
for(j in 1:num.y) {
if(subHazard) {
id <- which(sort.data$time>=y[j]|sort.data$index %w/o% c(0,event))
atRisk[id, j] <- G(y[j])/G(pmin(y[j], sort.data$time[id]))
} else {
id <- which(sort.data$time >= y[j])
atRisk[id, j] <- 1
}
}
return(list(y=y, d=d, delta=delta, atRisk=atRisk))
}
atRisk <- lapply(1:m, function(x) prepData(x)$atRisk)
delta <- lapply(1:m, function(x) prepData(x)$delta)
y <- lapply(1:m, function(x) prepData(x)$y)
d <- lapply(1:m, function(x) prepData(x)$d)



# Variables for saving the iteration history of the parameters
beta <- array(dim=c(p,m,MAX.ITER+1),dimnames=list(colnames(X),Type=1:m,NULL))
theta <- array(dim=c(min(m,s), min(m,s), MAX.ITER+1))
v <- array(dim=c(n, min(m,s), MAX.ITER+1),
dimnames=list(Group=1:n, Type=1:min(m,s), NULL))
# Initial values
if(!all(c("beta","theta","v")%in%names(inits))) stop("Invalid initial values")
if(min(m,s)!=ncol(as.matrix(inits$v))) stop("Invalid initial values")
beta[,,1] <- inits$beta
theta[,,1] <- inits$theta
v[,,1] <- inits$v
# Newton-Raphson Method
# -----------------------------------------------------------------------------
exp.eta <- lambda.0 <- Lambda.0 <- list()
diag.exp.eta <- diag.Lambda.0 <- lambda.0.sq <- list()
dexp.eta <- dlambda.0 <- dLambda.0 <- list()
ddiag.exp.eta <- ddiag.Lambda.0 <- dlambda.0.sq <- list()
d2exp.eta <- d2lambda.0 <- d2Lambda.0 <- list()
d2diag.exp.eta <- d2diag.Lambda.0 <- d2lambda.0.sq <- list()
ddiag.exp.eta<-ddiag.Lambda.0<-dlambda.0<-dlambda.0.sq<-replicate(s,list())
W <- replicate(m, list())
dW <- replicate(s, list())
d2W <- replicate(sum(1:s), list())
q.index = cbind(rep(1:s, times=s:1), unlist(sapply(1:s, function(q) q:s)))
colnames(q.index) = c("q1", "q2")
for(i in 1:MAX.ITER) {
SigmaInv = solve(theta[,,i])
SigmaInv.dSigma.SigmaInv = lapply(1:s, function(q)
solve(theta[,,i])%*%dSigma(q)%*%solve(theta[,,i]))
Q <- kronecker(SigmaInv, diag(n))
Q.prime<-lapply(SigmaInv.dSigma.SigmaInv,function(x) kronecker(-x,diag(n)))
Q.prime2 <- mapply(q1=q.index[,"q1"], q2=q.index[,"q2"], function(q1, q2)
kronecker(SigmaInv.dSigma.SigmaInv[[q2]]%*%dSigma(q1)%*%SigmaInv +
SigmaInv%*%dSigma(q1)%*%SigmaInv.dSigma.SigmaInv[[q2]],diag(n)),
SIMPLIFY=FALSE)



for(k in 1:m) {
exp.eta[[k]] <- as.vector(exp(X%*%beta[,k,i] + Z%*%v[,min(k,s),i]))
lambda.0[[k]] <- as.vector(d[[k]])/as.vector(t(atRisk[[k]])%*%
exp.eta[[k]])
Lambda.0[[k]] <- as.vector(atRisk[[k]]%*%diag(lambda.0[[k]])%*%
rep(1, ncol(atRisk[[k]])))
diag.exp.eta[[k]] <- diag(exp.eta[[k]])
diag.Lambda.0[[k]] <- diag(Lambda.0[[k]])
lambda.0.sq[[k]] <- diag(lambda.0[[k]]^2/d[[k]])
W[[k]] <- diag.exp.eta[[k]]%*%diag.Lambda.0[[k]] -
diag.exp.eta[[k]]%*%atRisk[[k]]%*%lambda.0.sq[[k]]%*%
t(atRisk[[k]])%*%diag.exp.eta[[k]]
}
# Gradient of profile h-like with respect to beta and v condtional on theta
# --------------------------------------------------------------------------
dhp.dbeta <- unlist(lapply(1:m, function(k) t(X)%*%
(delta[[k]]-diag.exp.eta[[k]]%*%Lambda.0[[k]])))
dhp.dv <- do.call(ifelse(frailty.cov=="none", "p.sum", "c"),
lapply(1:m, function(k) t(Z)%*%(delta[[k]]-
diag.exp.eta[[k]]%*%Lambda.0[[k]])))
G <- c(dhp.dbeta, dhp.dv + kronecker(-SigmaInv,diag(n))%*%as.vector(v[,,i]))
# Hessian of profile h-like with respect to beta and v condtional on theta
# --------------------------------------------------------------------------
H.XX <- block(lapply(W, function(w) t(X)%*%w%*%X))
H.XZ <- do.call(ifelse(frailty.cov=="none", "rbind", "block"),
lapply(W, function(w) t(X)%*%w%*%Z))
H.ZX <- do.call(ifelse(frailty.cov=="none", "cbind", "block"),
lapply(W, function(w) t(Z)%*%w%*%X))
H.ZZ <- do.call(ifelse(frailty.cov=="none", "p.sum", "block"),
lapply(W, function(w) t(Z)%*%w%*%Z))
H <- rbind(cbind(H.XX, H.XZ),
cbind(H.ZX, H.ZZ + Q))
H.inv <- solve(H)

################## Add ###############
H0 <- rbind(cbind(H.XX, H.XZ),
 cbind(H.ZX, H.ZZ))



# Update estimates of beta and v
# -------------------------------------------------------------------------
est <- c(beta[,,i], v[,,i]) + H.inv%*%G
beta[,,i+1] <- est[1:(m*p),1]
v[,,i+1] <- est[-(1:(m*p)),1]
# Gradient of the adjusted profile h-like with respect to theta given
# beta and v
# -------------------------------------------------------------------------
dv.dtheta <- lapply(1:s,function(q)
matrix(-solve(H.ZZ+Q)%*%Q.prime[[q]]%*%as.vector(v[,,i]),ncol=min(m,s)))
if(order>0)
for(q in 1:s) {
for(k in 1:m) {
dexp.eta[[k]] <- as.vector(exp.eta[[k]]*
(Z%*%dv.dtheta[[q]][,min(k,s)]))
dlambda.0[[q]][[k]] <- -(d[[k]]/as.vector((t(atRisk[[k]])%*%
exp.eta[[k]])^2))*as.vector(t(atRisk[[k]])%*%dexp.eta[[k]])
dLambda.0[[k]] <- as.vector(atRisk[[k]]%*%
diag(as.vector(dlambda.0[[q]][[k]]))%*%rep(1,ncol(atRisk[[k]])))
ddiag.exp.eta[[q]][[k]] <- diag(dexp.eta[[k]])
ddiag.Lambda.0[[q]][[k]] <- diag(dLambda.0[[k]])
dlambda.0.sq[[q]][[k]] <- diag(2*lambda.0[[k]]*
dlambda.0[[q]][[k]]/d[[k]])
diag.exp.eta.atRisk = diag.exp.eta[[k]]%*%atRisk[[k]]
dW[[q]][[k]] <- (ddiag.exp.eta[[q]][[k]]*diag.Lambda.0[[k]])+
(diag.exp.eta[[k]]*ddiag.Lambda.0[[q]][[k]])-
(ddiag.exp.eta[[q]][[k]]%*%atRisk[[k]]%*%lambda.0.sq[[k]]%*%
t(diag.exp.eta.atRisk))-
(diag.exp.eta.atRisk%*%dlambda.0.sq[[q]][[k]]%*%
t(diag.exp.eta.atRisk)) -
(diag.exp.eta.atRisk%*%lambda.0.sq[[k]]%*%
t(ddiag.exp.eta[[q]][[k]]%*%atRisk[[k]]))
}
}
if(order==0) dW<-replicate(s,list(replicate(m,list(matrix(0,nrow=N,ncol=N)))))




dH.XX <- lapply(dW,function(dW) block(lapply(dW,function(w) t(X)%*%w%*%X)))
dH.XZ <- lapply(dW, function(dW)
do.call(ifelse(frailty.cov=="none", "rbind", "block"),
lapply(dW, function(w) t(X)%*%w%*%Z)))
dH.ZX <- lapply(dW, function(dW)
do.call(ifelse(frailty.cov=="none", "cbind", "block"),
lapply(dW, function(w) t(Z)%*%w%*%X)))
dH.ZZ <- lapply(dW, function(dW)
do.call(ifelse(frailty.cov=="none", "p.sum", "block"),
lapply(dW, function(w) t(Z)%*%w%*%Z)))
dH.dtheta <- lapply(1:s, function(q)
rbind(cbind(dH.XX[[q]], dH.XZ[[q]]),
cbind(dH.ZX[[q]], dH.ZZ[[q]] + Q.prime[[q]])))
G.theta <- sapply(1:s, function(q)
sum(apply(as.matrix(v[,,i]), 1, function(v) -
.5*trace(SigmaInv%*%dSigma(q))+
.5*t(v)%*%SigmaInv.dSigma.SigmaInv[[q]]%*%v))-
.5*trace(H.inv%*%dH.dtheta[[q]]))
# Hessian of the adjusted profile h-like with respect to theta given
# beta and v
# --------------------------------------------------------------------
dv2.dtheta2 <- lapply(1:sum(1:s), function(q) {
q1 = q.index[q,"q1"]
q2 = q.index[q,"q2"]
matrix(-solve(H.ZZ+Q)%*%((dH.ZZ[[q1]] + Q.prime[[q1]])%*%
as.vector(dv.dtheta[[q2]]) + Q.prime[[q2]]%*%
as.vector(dv.dtheta[[q1]])+Q.prime2[[q]]%*%
as.vector(v[,,i])), ncol=min(m,s))
})
if(order>0)
for(q in 1:sum(1:s)) {
q1=q.index[q,"q1"]
q2=q.index[q,"q2"]
for(k in 1:m) {
d2exp.eta[[k]] <- as.vector((Z%*%dv.dtheta[[q1]][,min(k,s)])*
(Z%*%dv.dtheta[[q2]][,min(k,s)])*
exp.eta[[k]]+(Z%*%dv2.dtheta2[[q]][,min(k,s)])*exp.eta[[k]])



d2lambda.0[[k]] <- as.vector(-(lambda.0.sq[[k]]%*%t(atRisk[[k]])%*%
ddiag.exp.eta[[q2]][[k]]%*%
Z%*%dv.dtheta[[q1]][,min(k,s)] +
lambda.0.sq[[k]]%*%t(atRisk[[k]])%*%diag.exp.eta[[k]]%*%
Z%*%dv2.dtheta2[[q]][,min(k,s)] +
dlambda.0.sq[[q2]][[k]]%*%t(atRisk[[k]])%*%diag.exp.eta[[k]]%*%
Z%*%dv.dtheta[[q1]][,min(k,s)]))
d2Lambda.0[[k]] <- as.vector(atRisk[[k]]%*%
diag(as.vector(d2lambda.0[[k]]))%*%rep(1, ncol(atRisk[[k]])))
d2diag.exp.eta[[k]] <- diag(d2exp.eta[[k]])
d2diag.Lambda.0[[k]] <- diag(d2Lambda.0[[k]])
d2lambda.0.sq[[k]] <- diag((2*lambda.0[[k]]*dlambda.0[[q1]][[k]]*
dlambda.0[[q2]][[k]] + 2*lambda.0[[k]]*d2lambda.0[[k]])/d[[k]])
diag.exp.eta.atRisk = diag.exp.eta[[k]]%*%atRisk[[k]]
ddiag.exp.eta.q1.atRisk = ddiag.exp.eta[[q1]][[k]]%*%atRisk[[k]]
ddiag.exp.eta.q2.atRisk = ddiag.exp.eta[[q2]][[k]]%*%atRisk[[k]]
aa = d2diag.exp.eta[[k]]*diag.Lambda.0[[k]] +
(ddiag.exp.eta[[q1]][[k]]*ddiag.Lambda.0[[q2]][[k]] +
ddiag.exp.eta[[q2]][[k]]*ddiag.Lambda.0[[q1]][[k]]) +
diag.exp.eta[[k]]*d2diag.Lambda.0[[k]]
bb = (d2diag.exp.eta[[k]]%*%atRisk[[k]]%*%lambda.0.sq[[k]]%*%
t(diag.exp.eta.atRisk)) +
(ddiag.exp.eta.q1.atRisk%*%dlambda.0.sq[[q2]][[k]]%*%
t(diag.exp.eta.atRisk)) +
(ddiag.exp.eta.q1.atRisk%*%lambda.0.sq[[k]]%*%
t(ddiag.exp.eta.q2.atRisk))
cc = (ddiag.exp.eta.q2.atRisk%*%dlambda.0.sq[[q1]][[k]]%*%
t(diag.exp.eta.atRisk)) +
(diag.exp.eta.atRisk%*%d2lambda.0.sq[[k]]%*%
t(diag.exp.eta.atRisk)) +
(diag.exp.eta.atRisk%*%dlambda.0.sq[[q1]][[k]]%*%
t(ddiag.exp.eta.q2.atRisk))
dd = (ddiag.exp.eta.q2.atRisk%*%lambda.0.sq[[k]]%*%
t(ddiag.exp.eta.q1.atRisk)) +
(diag.exp.eta.atRisk%*%dlambda.0.sq[[q2]][[k]]%*%
t(ddiag.exp.eta.q1.atRisk)) +
(diag.exp.eta.atRisk%*%lambda.0.sq[[k]]%*%
t(d2diag.exp.eta[[k]]%*%atRisk[[k]]))
d2W[[q]][[k]] = aa-(bb+cc+dd)
}
}
if(order==0) d2W <- replicate(sum(1:s),
list(replicate(m, list(matrix(0, nrow=N, ncol=N)))))



d2H.XX <- lapply(d2W, function(d2W)
block(lapply(d2W, function(w) t(X)%*%w%*%X)))
d2H.XZ <- lapply(d2W, function(d2W)
do.call(ifelse(frailty.cov=="none", "rbind", "block"),
lapply(d2W, function(w) t(X)%*%w%*%Z)))
d2H.ZX <- lapply(d2W, function(d2W)
do.call(ifelse(frailty.cov=="none", "cbind", "block"),
lapply(d2W, function(w) t(Z)%*%w%*%X)))
d2H.ZZ <- lapply(d2W, function(d2W)
do.call(ifelse(frailty.cov=="none", "p.sum", "block"),
lapply(d2W, function(w) t(Z)%*%w%*%Z)))
d2H.dtheta2 <- lapply(1:sum(1:s), function(q) {
rbind(cbind(d2H.XX[[q]], d2H.XZ[[q]]),
cbind(d2H.ZX[[q]], d2H.ZZ[[q]] + Q.prime2[[q]]))})
H.theta.temp <- sapply(1:sum(1:s), function(q)
{q1=q.index[q,"q1"]; q2=q.index[q,"q2"]
sum(apply(as.matrix(v[,,i]), 1, function(v) {
.5*trace(SigmaInv%*%dSigma(q2)%*%SigmaInv%*%dSigma(q1)) +
.5*t(v)%*%(SigmaInv%*%dSigma(q2)%*%SigmaInv.dSigma.SigmaInv[[q1]])%*%v +
.5*t(v)%*%(SigmaInv%*%dSigma(q1)%*%SigmaInv.dSigma.SigmaInv[[q2]])%*%v}))-
c(v[,,i])%*%Q.prime[[q1]]%*%c(dv.dtheta[[q2]]) +
.5*trace(-H.inv%*%dH.dtheta[[q1]]%*%H.inv%*%dH.dtheta[[q2]]+
H.inv%*%d2H.dtheta2[[q]])})
H.theta <- matrix(nrow=s, ncol=s)
H.theta[lower.tri(H.theta, diag=TRUE)] <- H.theta.temp
H.theta[upper.tri(H.theta)] <- H.theta[lower.tri(H.theta)]




# Update estimates of theta
# -------------------------------------------------------------------------
if(order==2 & s!=1) stop("Invalid order argument")
if(order==0 | order==1) {
if(frailty.cov=="independent") {
theta[,,i+1] <- diag(as.vector(diag(theta[,,i]) +
solve(H.theta)%*%G.theta))
} else if(frailty.cov=="exchangeable") {
temp <- c(diag(theta[,,i]),unique(theta[,,i][lower.tri(theta[,,i])]))+
solve(H.theta)%*%G.theta
theta[,,i+1] <- diag(temp[1:(s-1)])
theta[,,i+1][lower.tri(theta[,,i+1])|upper.tri(theta[,,i+1])]<-temp[s]
} else if(frailty.cov=="unstructured") {
temp <- theta[,,i][lower.tri(theta[,,i], diag=TRUE)]+
solve(H.theta)%*%G.theta
theta[,,i+1][lower.tri(theta[,,i], diag=TRUE)] <- temp
theta[,,i+1][upper.tri(theta[,,i+1], diag=TRUE)] <-
theta[,,i+1][lower.tri(theta[,,i+1], diag=TRUE)]
} else if(frailty.cov=="none") {
theta[,,i+1] <- theta[,,i] + solve(H.theta)%*%G.theta
}
} else {
delta.i.plus <- lapply(delta, function(x)
tapply(x, sort.data$cluster, sum))
mu.i.plus <- lapply(1:m, function(k)
tapply(Lambda.0[[k]]*exp(X%*%beta[,k,i]), sort.data$cluster, sum))
dh.dv <- function(v, j) {
Reduce("+",delta.i.plus)[j]-exp(v)*
Reduce("+",mu.i.plus)[j]-v/theta[,,i]}
v.tilda <- sapply(1:n, function(x)
uniroot(dh.dv, c(-50,50), j=x, tol=1E-3)$root)
A <- exp(v.tilda)*Reduce("+", mu.i.plus)
trace.dS <-sum(6*A*theta[,,i]^(-2)*(A+1/theta[,,i])^(-3)-15*
A^2*theta[,,i]^(-2)*(A+1/theta[,,i])^(-4))
trace.d2S<-sum(12*A*theta[,,i]^(-3)*(A+1/theta[,,i])^(-3)-18*
A*theta[,,i]^(-4)*(A+1/theta[,,i])^(-4) -
30*A^2*theta[,,i]^(-3)*(A+1/theta[,,i])^(-4)+60*
A^2*theta[,,i]^(-4)*(A+1/theta[,,i])^(-5))
G.theta.order2 <- G.theta - trace.dS/24
H.theta.order2 <- H.theta - trace.d2S/24
theta[,,i+1] <- theta[,,i] + G.theta.order2/H.theta.order2
}




# Convergence
# -------------------------------------------------------------------------
converge=0
if(max(abs(c(beta[,,i+1], theta[,,i+1])-
c(beta[,,i], theta[,,i])), na.rm=TRUE) < TOL){
converge=1
# break
}
}
# Correction: Profile h-likelihood and Adjusted Profile h-likelihood
# -------------------------------------------------------------------------
h.1<- sum(sapply(1:m, function(k) sum(delta[[k]]*log(exp.eta[[k]])))) -
   sum(sapply(1:m, function(k) sum(d[[k]]*log(as.vector(t(atRisk[[k]])%*%exp.eta[[k]])))))

 h.p <- sum(sapply(1:m, function(k) sum(delta[[k]]*log(exp.eta[[k]])))) - 
		sum(sapply(1:m, function(k) sum(d[[k]]*log(as.vector(t(atRisk[[k]])%*%exp.eta[[k]]))))) +
		.5*log(det(Q/2/pi))-.5*t(c(v[,,i+1]))%*%Q%*%c(v[,,i+1])


#h.2<- sum(apply(as.matrix(v[,,i]), 1, function(v)
#.5*log(det(SigmaInv/2/pi))-.5*t(v)%*%SigmaInv%*%v))
#h.p <- h.1 + h.2
#A.h.p <- h.p - .5*log(det(H))+(ncol(X)+ncol(Z))/2*log(2*pi)

A.h.p <- h.p - .5*log(det(H))+ ncol(H)/2*log(2*pi)

#A.h.p <- h.p - .5*log(det(H))+(p*m+n)/2*log(2*pi)

rAIC <- -2*A.h.p + 2*s # Akaike information criterion

cAIC<- -2*h.1 +2*trace(solve(H)%*%H0)

# Profile h-likelihood and Adjusted Profile h-likelihood
# -------------------------------------------------------------------------
#h.1<- sum(sapply(1:m, function(k) sum(delta[[k]]*log(exp.eta[[k]])))) +
#- sum(sapply(1:m, function(k)
#sum(d[[k]]*log(as.vector(t(atRisk[[k]])%*%exp.eta[[k]])))))

#h.p <- sum(sapply(1:m, function(k) sum(delta[[k]]*log(exp.eta[[k]])))) +
#- sum(sapply(1:m, function(k)
#sum(d[[k]]*log(as.vector(t(atRisk[[k]])%*%exp.eta[[k]]))))) +
#sum(apply(as.matrix(v[,,i]), 1, function(v)
#.5*log(det(SigmaInv/2/pi))-.5*t(v)%*%SigmaInv%*%v))

#A.h.p <- h.p - .5*log(det(H))+(p+n)/2*log(2*pi)
#AIC <- -2*A.h.p + 2*s # Akaike information criterion


# Confidence intervals and hypothesis tests
# -------------------------------------------------------------------------
# Estimates, standard errors, and CI for regression coefficients
beta.SE <- sqrt(diag(as.matrix(H.inv[1:(m*p), 1:(m*p)])))
lower.beta <- c(beta[,,i]) - qnorm(1-alpha/2)*beta.SE
upper.beta <- c(beta[,,i]) + qnorm(1-alpha/2)*beta.SE
beta.CI <- data.frame(rep(1:m,each=p), rep(colnames(X),times=m),
c(beta[,,i]), beta.SE, lower.beta, upper.beta)
colnames(beta.CI) <- c("Type", "Effect", "Estimate", "SE",
paste(alpha/2*100, "%", sep=""), paste((1-alpha/2)*100, "%", sep=""))
# Variance component estimates
Var.Comp <- paste("Sigma.", unlist(sapply(1:min(m,s),
function(q) q:min(m,s))), rep(1:min(m,s), times=min(m,s):1), sep="")
theta.Est <- data.frame(Var.Comp,
Estimate=theta[,,i][lower.tri(theta[,,i], diag=TRUE)])



# Test H0: theta=0 for cause-specific and subHazard models when univariate
# normal distribution
if(subHazard==TRUE) {
loglik.noFrailty <- crr(sort.data$time, sort.data$index, X)$loglik
} else {
loglik.noFrailty <- sum(sapply(1:m,function(k)
coxph(Surv(time,index==k)~X,data=sort.data)$loglik[2]))
}
if(frailty.cov=="none") {
psi.hat <- -2*(loglik.noFrailty - A.h.p)
theta.pvalue <- ifelse(psi.hat<=0,1,0.5*pchisq(psi.hat,1,lower.tail=FALSE))
} else theta.pvalue=NULL
# Gather and output results
# -------------------------------------------------------------------------

######################### Correction ##########################################
v.out <- Z%*%v[,,i]; rownames(v.out)<-sort.data$cluster;
out <- list(beta=beta.CI, theta=theta.Est, v=v.out, theta.pvalue=theta.pvalue,
lambda.0=lapply(1:m, function(k) matrix(lambda.0[[k]],
dimnames=list(Time=y[[k]], Type=k))),
time=y, h1=h.1, hp=h.p, Ahp=A.h.p, loglik.noFrailty=loglik.noFrailty, rAIC=rAIC, cAIC=cAIC, H=H,H0=H0,
#time=y, h1=h.1, h2=h.2, hp=h.p, Ahp=A.h.p, loglik.noFrailty=loglik.noFrailty, rAIC=rAIC, cAIC=cAIC, H=H,H0=H0,
converge=converge, iterations=i, call=match.call())
class(out) <- "hlike"
print(out$beta); cat("\n"); print(out$theta);
cat("\n",ifelse(out$converge==1, "Successfully Converged",
"Failed to Converge"),"\n")
return(out)
}
