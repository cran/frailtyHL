jmfit2 <-
function(jm1,jm2,data_conti,data_surv, Maxiter=200) {

mc <- match.call()
res1<-design_hglm(jm1[[3]],DataMain=data_conti)
res2<-design_frailtyHL(jm2[[3]],data=data_surv)

rho<- 0
phi<- 1
alpha<- 0.1

## x matrix
xx<-dbind(res1[[2]],res2[[2]])
## z matrix
nrow(res1[[3]])
zz<-rbind(res1[[3]],rho*res2[[4]][[1]])
pp<-ncol(xx)
qq<-ncol(zz)
namesX1<-res1[[8]]
namesX2<-res2[[15]]

x<-res2[[2]]
beta_h0<-res2[[9]]
z<-res2[[4]][[1]]
v_h0<-res2[[11]]
di<-res2[[8]]
Mi<-res2[[5]]
idx2<-res2[[6]]

convergence<-1
i<-0

iteration<-0
#Maxiter=5

while (convergence>0.0001  && iteration<=Maxiter) {
 # while (convergence>0.001  && iteration<=Maxiter) {
iteration<-iteration+1
#while (convergence>0.001) {

## weight matrix
zz<-rbind(res1[[3]],rho*res2[[4]][[1]])
weight1<-diag(rep(1,res1[[5]])/phi)
muh <- x %*% beta_h0 + rho*z %*% v_h0
expeta <- exp(muh)
cla0 <- di/(crossprod(Mi,expeta))
Ai <- diag(cla0[, 1])
Wi <- diag(expeta[, 1])
done <- matrix(1, idx2, 1)
clam0 <- Mi %*% (diag(Ai)* done)
Bi <- diag(clam0[, 1])
temp4 <- cla0^2/di
Adi <- diag(temp4[, 1])


weight2<-(diag(Wi) * Bi) - (diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))
weight3<-diag(rep(1/alpha,qq))

weight<-dbind(dbind(weight1,weight2),weight3)
weight2_star<-diag(rep(1,nrow(weight2)))
weight_star<-dbind(dbind(weight1,weight2_star),weight3)

####
oo<-matrix(0,qq,pp)
ii<-diag(rep(1,qq))
TT<-rbind(cbind(xx,zz),cbind(oo,ii))

### adjusted dependent
adj_zz1<-res1[[1]]
adj_zz2<-weight2%*%muh+(res2[[3]]-clam0*expeta)
adj_zz3<-matrix(rep(0,qq),qq,1)
adj_zz<-rbind(adj_zz1,adj_zz2,adj_zz3)
hinv<-solve(t(TT)%*%weight%*%TT)

beta_vv<-hinv%*%t(TT)%*%weight_star%*%adj_zz  ## Estimation of (beta,v)

temp1<-res1[[6]]+1
temp2<-res1[[6]]+res2[[12]]
beta_h0<-matrix(beta_vv[temp1:temp2,1],res2[[12]],1)
temp1<-1+pp
temp2<-res1[[7]]+pp
v_h0<-matrix(beta_vv[temp1:temp2,1],qq,1)

old_phi<-phi
old_alpha<-alpha

n1<-res1[[5]]
p1<-res1[[6]]
mu1<-res1[[2]]%*%beta_vv[1:p1]+res1[[3]]%*%v_h0
weight_star<- dbind(dbind(-1*diag(rep(1,n1)/phi^2),0*weight2),0*weight3)
gamma0<- phi*sum(hinv*(t(TT)%*%weight_star%*%TT))
phi<-t(res1[[1]]-mu1)%*%(res1[[1]]-mu1)/(n1-(p1+qq-gamma0))
phi<-t(res1[[1]]-mu1)%*%(res1[[1]]-mu1)/(n1+gamma0) ## Estimation of phi=var(error)

weight_star2<- dbind(dbind(0*weight1,0*weight2),-1*diag(rep(1/alpha^2,qq)))
gamma1<- -alpha*sum(hinv*(t(TT)%*%weight_star2%*%TT))
alpha<-t(v_h0)%*%(v_h0)/(qq-gamma1)   ## Estimation of alpha=var(v)

zz_rho<-rbind(0*res1[[3]],res2[[4]][[1]])
TT1<-rbind(cbind(0*xx,zz_rho),cbind(oo,0*ii))
Hda<- 2*t(TT1)%*%weight%*%TT


#dexpeta=expeta*(z%*%v_h0)
#dm <- (crossprod(Mi,dexpeta))/(crossprod(Mi,expeta))
#cla0_rho <- -cla0*dm

cla0_rho <- -cla0

Ai_rho <- diag(cla0_rho[,1])
clam0_rho <- Mi %*% (diag(Ai_rho)* done)
Bi_rho <- diag(clam0_rho[, 1])

Adi_rho<- (2*cla0/di)*cla0_rho
#Adi_rho<- 2*cla0/di*cla0_rho

Adi_rho<-diag(Adi_rho[,1])


## weight2<-(diag(Wi) * Bi) - (diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))

weight2_rho<-weight2 * diag(z%*%v_h0)

weight2_rho<-weight2_rho+(diag(Wi) * Bi_rho) * diag(z%*%v_h0)

weight2_rho<-weight2_rho-(diag(Wi) * Mi) %*% (diag(Adi_rho)*( t(diag(Wi)*Mi))) * diag(z%*%v_h0)

weight2_rho<-weight2_rho-(diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))*diag(z%*%v_h0)

weight_rho<-dbind(dbind(0*weight1,weight2_rho),0*weight3)

dhdrho <- sum((res2[[3]]-clam0*expeta) * (z %*% v_h0))
dhdrho <- dhdrho-0.5*sum(hinv*Hda)-0.5*sum(hinv*(t(TT)%*%weight_rho%*%TT))

d2hdrho2 <- t(z %*% v_h0) %*% (diag(Wi)* Bi) %*% (z %*% v_h0)

old_rho<-rho
rho<-as.vector(rho+dhdrho/d2hdrho2)

convergence<-abs(old_phi-phi)+abs(old_alpha-alpha)+abs(old_rho-rho)
i<-i+1
####################################################################

beta_h<-beta_vv[1:pp]
se_betah<-sqrt(diag(hinv))
se_beh<-se_betah[1:pp]
t_value=beta_h/se_beh
p_value=2*(1-pnorm(abs(t_value)))

rho_h<-rho
phi_h<-as.vector(phi)
alpha_h<-as.vector(alpha)




#print("v est : ")
#print(v_h0)

#print("beta est : ")
#print(beta_h)
#print("SE of beta_h : ")
#print(se_beh)
#print("p-value")
#print(p_value)
#print("phi : ")
#print(phi)
#print("alpha : ")
#print(alpha)
#print("rho : ")
#print(rho)
#print("convergence : ")
#print(convergence)
}

#print("iteration : ")
#print(iteration)
print("iterations : ")
print(i)
print("convergence : ")
print(convergence)


F.Est<-cbind(beta_h,se_beh, t_value, p_value)  # Est. of fixed effects
D.Est<-cbind(phi_h,alpha_h,rho_h)     # Est. of dispersion paras
V.Est<-as.vector(v_h0)

#print(F.Est)
#print(D.Est)
#print(V.Est)

#res<-list(F.Est=F.Est,D.Est=D.Est)
res<-list(F.Est=F.Est,D.Est=D.Est,V.Est=V.Est)

print("Estimates for fixed effects")
print_mat1<-res$F.Est
colnames(print_mat1)<-c("Estimate","Std. Error", "t_value", "p_value")
rownames(print_mat1) <- c(namesX1,namesX2)
print(round(print_mat1,digits=5))
#print(print_mat1,4)

print("Estimates for dispersion parameters")
print_mat2<-res$D.Est
colnames(print_mat2)<-c("phi_h", "alpha_h", "gamma_h" )
print(round(print_mat2,digits=5))
#print(print_mat2,4)
return(res)
}
