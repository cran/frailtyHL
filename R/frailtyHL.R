frailtyHL <-
function(formulaMain,censor,DataMain,RandDist="Normal",mord=0,dord=1,Maxiter=200,convergence=10^-6){
    require(Matrix)
    require(numDeriv)
    mc <- match.call()
    fr <- HGLMFrames(mc, formulaMain, contrasts)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formulaMain, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(fr$Y, length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <-matrix(0,n,1)
    del[,1] <- censor
    SS <- FL$Subject
    res1<-FrailtyMakeData(y,x,del,z)
    y<-res1[1][[1]]
    x<-res1[2][[1]]
    del<-res1[3][[1]]
    z<-res1[4][[1]]
    Mi<-res1[5][[1]]
    idx2<-res1[6][[1]]
    t2<-res1[7][[1]]
    di<-res1[8][[1]]
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    alpha_h <- rep(0, nrand)
    for (i in 1:nrand) alpha_h[i] <- 0.1
    Max_iter<-Maxiter
    err<-1
    for ( i in 1:Max_iter) {
        if (err>=convergence) {
        if (RandDist=="Normal") res2<-PNfrailtyHL(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord)
        if (RandDist=="Gamma") res2<-PGfrailtyHL(x,z,y, del,Mi,idx2,t2, di, beta_h,v_h, alpha_h,mord,dord)
        alpha_h<-res2[13][[1]]
        alpha_h1<-res2[14][[1]]
        temp4<-sum(abs(alpha_h-alpha_h1))
        err<-temp4
        beta_h<-res2[11][[1]]
        v_h<-res2[12][[1]]
        alpha_h<-alpha_h1
        se_beta<-res2[20][[1]]
        print_i<-i
        print_err<-err
        names(print_i) <- "iteration : "
        print(print_i)
        names(print_err) <- "convergence : "
        print(print_err)
        }
    }
    if (err<convergence) print("converged") 
    if (err>convergence) print("did not converge")
###############################################################
############# print estimates ###########################
###############################################################
    result<-list(0)
    names(result)[1]<-"Model"
    if (RandDist=="Gamma") {
         print("Results from the gamma frailty model")
         result$Model<-"gamma frailty model"
    }
    if (RandDist=="Normal") {
         print("Results from the log-normal frailty model")
         result$Model<-"log-normal frailty model"
    }
    result$formulaMain<-formulaMain
    print(formulaMain)
    if (mord==0 && dord==1) {
         print("Method : HL(0,1)")   
         result$Method<-"HL(0,1)"
    }
    if (mord==0 && dord==2) {
         print("Method : HL(0,2)")     
         result$Method<-"HL(0,2)"
    }
    if (mord==1 && dord==1) {
         print("Method : HL(1,1)")   
         result$Method<-"HL(1,1)"
    }
    if (mord==1 && dord==2) {
         print("Method : HL(1,2)")
         result$Method<-"HL(1,2)"
    }
    print("Estimates from the mean model")
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
    result$FixCoef<-beta_coeff
    print(beta_coeff,4)
###############################################################
############# se for lambda         ###########################
###############################################################
    if (RandDist=="Normal") res3<-PNFrailty_SE.h(res2,nrand,q,qcum,dord)
    if (RandDist=="Gamma") res3<-PGFrailty_SE.h(res2,nrand,q,qcum,dord)
    print("Estimates from the dispersion model")
    se_alpha_h<-res3[1][[1]]
    hlike<--2*res3[2][[1]]
    p1<--2*res3[3][[1]]
    p2<--2*res3[4][[1]]
    p3<--2*res3[5][[1]]
    z_lam<-alpha_h/se_alpha_h
    lam_coeff<-cbind(alpha_h,se_alpha_h)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    print(lam_coeff,4)
    result$RandCoef<-lam_coeff
###############################################################
############# -2*Likelihoods         ##########################
###############################################################
    if (mord==0 && dord==1) like_value<-cbind(hlike,p1)
    if (mord==0 && dord==1) colnames(like_value) <- c("-2*hp","-2*p_b,v(hp)")
    if (mord==0 && dord==2) like_value<-cbind(hlike,p3)
    if (mord==0 && dord==2) colnames(like_value) <- c("-2*hp","-2*s_b,v(hp)")
    if (mord==1 && dord==1) like_value<-cbind(hlike,p2,p1)
    if (mord==1 && dord==1) colnames(like_value) <- c("-2*hp","-2*p_v(hp)","-2*p_b,v(hp)")
    if (mord==1 && dord==2) like_value<-cbind(hlike,p2,p3)
    if (mord==1 && dord==2) colnames(like_value) <- c("-2*hp","-2*p_v(hp)","-2*s_b,v(hp)")
    result$likelihood<-like_value
    result$iter<-print_i
    if (print_err<convergence) result$convergence<-"converged"
    if (print_err>convergence) result$convergence<-"did not converge"
    names(result$convergence) <- "convergence : "
    print(like_value,5)
    res4<-list(res2,res3)
    return(result)   
}

