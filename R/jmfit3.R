jmfit3 <-
function(jm1,jm2,jm3,data_conti,data_surv,Maxiter=200) {
mc<-match.call()
res1<-design_hglm(jm1[[3]],DataMain=data_conti)
res2<-design_frailtyHL(jm2[[3]],data=data_surv)
res3<-design_frailtyHL(jm3[[3]],data=data_surv)

 rho<- 1
 rho2<- 1

 #rho<- -10
 #rho2<- -10

phi<- 1
alpha<- 0.1
namesX1<-res1[[8]]
namesX2<-res2[[15]]
namesX3<-res3[[15]]


## x matrix
xx<-dbind(dbind(res1[[2]],res2[[2]]),res3[[2]])
## z matrix
nrow(res1[[3]])
zz<-rbind(res1[[3]],rho*res2[[4]][[1]], rho2*res3[[4]][[1]])
pp<-ncol(xx)
qq<-ncol(zz)

x<-res2[[2]]
beta_h0<-res2[[9]]
z<-res2[[4]][[1]]
v_h0<-res2[[11]]
di<-res2[[8]]
Mi<-res2[[5]]
idx2<-res2[[6]]

x2<-res2[[2]]
beta_h02<-res3[[9]]
z2<-res3[[4]][[1]]
v_h02<-res3[[11]]
di2<-res3[[8]]
Mi2<-res3[[5]]
idx2_2<-res3[[6]]

convergence<-1
i<-0
iteration<-0

#Maxiter=200

while (convergence>0.0001  && iteration<=Maxiter) {
iteration<-iteration+1


## weight matrix
#zz<-rbind(res1[[3]],rho*res2[[4]][[1]])
zz<-rbind(res1[[3]],rho*res2[[4]][[1]], rho2*res3[[4]][[1]])
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
weight2<-(diag(Wi) * Bi) - (diag(Wi)* Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))


muh2 <- x2 %*% beta_h02 + rho2*z2 %*% v_h0

#muh2 <- x2 %*% beta_h02 + rho2*z2 %*% v_h02
expeta2 <- exp(muh2)
cla02 <- di2/(crossprod(Mi2,expeta2))
Ai2 <- diag(cla02[, 1])
Wi2 <- diag(expeta2[, 1])
done2 <- matrix(1, idx2_2, 1)
clam02 <- Mi2 %*% (diag(Ai2)* done2)
Bi2 <- diag(clam02[, 1])
temp4_2 <- cla02^2/di2
Adi2 <- diag(temp4_2[, 1])
weight2_2<-(diag(Wi2)*Bi2) - (diag(Wi2)* Mi2) %*% (diag(Adi2)*( t(diag(Wi2)*Mi2)))

weight3<-diag(rep(1/alpha,qq))  #### U

 weight<-dbind(dbind(dbind(weight1,weight2), weight2_2),weight3)
 weight2_star<-diag(rep(1,nrow(weight2)))
 weight3_star<-diag(rep(1,nrow(weight2_2)))
 weight_star<-dbind(dbind(dbind(weight1,weight2_star), weight3_star),weight3)

####
oo<-matrix(0,qq,pp)
ii<-diag(rep(1,qq))
TT<-rbind(cbind(xx,zz),cbind(oo,ii))

### adjusted dependent   
adj_zz1<-res1[[1]]
adj_zz2<-weight2%*%muh+(res2[[3]]-clam0*expeta)
adj_zz2_2<-weight2_2%*%muh2+(res3[[3]]-clam02*expeta2)

adj_zz3<-matrix(rep(0,qq),qq,1)
adj_zz<-rbind(adj_zz1,adj_zz2,adj_zz2_2, adj_zz3)
hinv<-solve(t(TT)%*%weight%*%TT)

beta_vv<-hinv%*%t(TT)%*%weight_star%*%adj_zz

beta_h<-beta_vv[1:pp]

temp1<-res1[[6]]+1
temp2<-res1[[6]]+res2[[12]]
beta_h0<-matrix(beta_vv[temp1:temp2,1],res2[[12]],1)

temp1_2<-res1[[6]]+res1[[6]]+1
temp2_2<-res1[[6]]+res2[[12]]+res3[[12]]
beta_h02<-matrix(beta_vv[temp1_2:temp2_2,1],res3[[12]],1)

temp1<-1+pp
temp2<-res1[[7]]+pp
v_h0<-matrix(beta_vv[temp1:temp2,1],qq,1)

######################### alpha  ########################################

old_phi<-phi
old_alpha<-alpha

n1<-res1[[5]]
p1<-res1[[6]]
mu1<-res1[[2]]%*%beta_vv[1:p1]+res1[[3]]%*%v_h0

weight_star<- dbind(dbind(dbind(-1*diag(rep(1,n1)/phi^2),0*weight2), 0*weight2_2),0*weight3)

gamma0<- phi*sum(hinv*(t(TT)%*%weight_star%*%TT))
#phi<-t(res1[[1]]-mu1)%*%(res1[[1]]-mu1)/(n1-(p1+qq-gamma0))
phi<-t(res1[[1]]-mu1)%*%(res1[[1]]-mu1)/(n1+gamma0)

weight_star2<- dbind(dbind(dbind(0*weight1,0*weight2), 0*weight2_2),-1*diag(rep(1/alpha^2,qq)))
gamma1<- -alpha*sum(hinv*(t(TT)%*%weight_star2%*%TT))
alpha<-t(v_h0)%*%(v_h0)/(qq-gamma1)


################################### rho ##################################################
 #zz<-rbind(res1[[3]],rho*res2[[4]][[1]], rho2*res3[[4]][[1]])

zz_rho<-rbind(0*res1[[3]],res2[[4]][[1]], 0*res3[[4]][[1]])
TT1<-rbind(cbind(0*xx,zz_rho),cbind(oo,0*ii))
Hda<- 2*t(TT1)%*%weight%*%TT

#dexpeta=expeta*(z%*%v_h0)
#dm <- (crossprod(Mi,dexpeta))/(crossprod(Mi,expeta))
#cla0_rho <- -cla0*dm

cla0_rho <- -cla0
Ai_rho <- diag(cla0_rho[,1])
clam0_rho <- Mi %*%(diag(Ai_rho)*done)
Bi_rho <- diag(clam0_rho[, 1])

Adi_rho<- (2*cla0/di)*cla0_rho
Adi_rho<-diag(Adi_rho[,1])

weight2_rho<-weight2 * diag(z%*%v_h0)

weight2_rho<-weight2_rho+(diag(Wi) * Bi_rho) * diag(z%*%v_h0)

weight2_rho<-weight2_rho-(diag(Wi) * Mi) %*% (diag(Adi_rho)*( t(diag(Wi)*Mi))) * diag(z%*%v_h0)

weight2_rho<-weight2_rho-(diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))*diag(z%*%v_h0)

weight_rho<-dbind(dbind(dbind(0*weight1,weight2_rho), 0*weight2_2),0*weight3)

dhdrho <- sum((res2[[3]]-clam0*expeta) * (z %*% v_h0))
dhdrho <- dhdrho-0.5*sum(hinv*Hda)-0.5*sum(hinv*(t(TT)%*%weight_rho%*%TT))

d2hdrho2 <- t(z %*% v_h0) %*% (diag(Wi)* Bi) %*% (z %*% v_h0)

old_rho<-rho
 rho<-as.vector(rho+dhdrho/d2hdrho2)



################################### rho2 ##################################################
 #zz<-rbind(res1[[3]],rho*res2[[4]][[1]], rho2*res3[[4]][[1]])

zz_rho2<-rbind(0*res1[[3]],0*res2[[4]][[1]], res3[[4]][[1]])
TT1_2<-rbind(cbind(0*xx,zz_rho2),cbind(oo,0*ii))
Hda2<- 2*t(TT1_2)%*%weight%*%TT

cla0_rho2 <- -cla02
Ai_rho2 <- diag(cla0_rho2[,1])
clam0_rho2 <- Mi2 %*% (diag(Ai_rho2)* done2)
Bi_rho2 <- diag(clam0_rho2[, 1])

Adi_rho2<- (2*cla02/di2)*cla0_rho2
Adi_rho2<-diag(Adi_rho2[,1])

weight2_rho2<-weight2_2 * diag(z2%*%v_h0)

weight2_rho2<-weight2_rho2+(diag(Wi2) * Bi_rho2) * diag(z2%*%v_h0)

weight2_rho2<-weight2_rho2-(diag(Wi2) * Mi2) %*% (diag(Adi_rho2)*( t(diag(Wi2)*Mi2))) * diag(z2%*%v_h0)

weight2_rho2<-weight2_rho2-(diag(Wi2) * Mi2) %*% (diag(Adi2)*( t(diag(Wi2)*Mi2)))*diag(z2%*%v_h0)

weight_rho2<-dbind(dbind(dbind(0*weight1, 0*weight2), weight2_rho2),0*weight3)

dhdrho_2 <- sum((res3[[3]]-clam02*expeta2) * (z2 %*% v_h0))
dhdrho_2 <- dhdrho_2-0.5*sum(hinv*Hda2)-0.5*sum(hinv*(t(TT)%*%weight_rho2%*%TT))

d2hdrho_2 <- t(z2 %*% v_h0) %*% (diag(Wi2)* Bi2) %*% (z2 %*% v_h0)

old_rho2<-rho2
 rho2<-as.vector(rho2+dhdrho_2/d2hdrho_2)

convergence<-abs(old_phi-phi)+abs(old_alpha-alpha)+abs(old_rho-rho)+abs(old_rho2-rho2)

i<-i+1
####################################################################

beta_h<-beta_vv[1:pp]
se_betah<-sqrt(diag(hinv))
se_beh<-se_betah[1:pp]
t_value=beta_h/se_beh
p_value=2*(1-pnorm(abs(t_value)))

rho_h<-rho
rho2_h<-rho2
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
#print("rho2 : ")
#print(rho2)
#print("convergence : ")
#print(convergence)
}


print("iterations : ")
print(i)

print("convergence : ")
print(convergence)

rho1_h<-rho_h

F.Est<-cbind(beta_h,se_beh, t_value, p_value)  # Est. of fixed effects
D.Est<-cbind(phi_h,alpha_h,rho1_h,rho2_h)      # Est. of dispersion paras
V.Est<-as.vector(v_h0)

#print(F.Est)
#print(D.Est)
#print(V.Est)

#res<-list(F.Est=F.Est,D.Est=D.Est)
res<-list(F.Est=F.Est,D.Est=D.Est,V.Est=V.Est)

print("Estimates for fixed effects")
print_mat3<-res$F.Est
colnames(print_mat3)<-c("Estimate","Std. Error", "t_value", "p_value")
rownames(print_mat3) <- c(namesX1,namesX2,namesX3)
print(round(print_mat3,digits=5))
#print(print_mat3,4)

print("Estimates for dispersion parameters")
print_mat4<-res$D.Est
colnames(print_mat4)<-c("phi_h", "alpha_h", "gamma1_h",  "gamma2_h")
print(round(print_mat4,digits=5))
#print(print_mat4,4)

return(res)

}
