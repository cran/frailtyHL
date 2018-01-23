AR1frailtyHL <-
function (x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, alpha_h0, 
    mord, dord, varfixed = FALSE,varnonneg,IArho,rho_0,qrep) 
{
    n <- nrow(x)
    p <- ncol(x)
    z <- diag(1,n,n)
    nrand <- 1
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z)[2]
    qcum <- cumsum(c(0, q))
    beta_h <- beta_h0
    v_h <- v_h0
    rho <- rho_0
    for (i in 1:nrand) {
        if (alpha_h0[i] < 1e-06) 
            alpha_h0[i] <- 1e-06
    }
    alpha_h <- alpha_h0
    zz <- z
    z <- matrix(0, n, qcum[nrand + 1])
    z[1:n, 1:qcum[nrand + 1]] <- zz[1:n, 1:qcum[nrand + 1]]
    muh <- x %*% beta_h0 + z %*% v_h0
    expeta <- exp(muh)
    cla0 <- di/(crossprod(Mi,expeta))
    Wi <- diag(expeta[, 1])
    Ai <- diag(cla0[, 1])
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    clam0 <- Mi %*% Ai %*% done
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    dft1 <- crossprod(x, del - Wi %*% clam0)
    if (mord == 1) {
        U <- iD
        Bi <- diag(clam0[, 1])
        temp4 <- cla0^2/di
        As <- diag(temp4[, 1])
        mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% As %*% (t(Mi) %*% Wi)
        dinv <- solve(t(z) %*% mat %*% z + U)
        mu0 <- exp(x %*% beta_h0 + z %*% v_h0)
        mu <- exp(x %*% beta_h0 + z %*% v_h0) * clam0
        dcla0be <- -As %*% (t(Mi) %*% Wi)
        dcla0b <- matrix(0, idx2, p)
        dv_db <- matrix(0, qcum[nrand + 1], p)
        xz <- matrix(0, n, p)
        dmu0 <- matrix(0, n, p)
        dw_db1 <- matrix(0, n, p)
        xk <- matrix(0, n, 1)
        ad1 <- matrix(0, p, 1)
        for (k in 1:p) {
            xk[, 1] <- x[, k]
            dv_db[, k] <- -dinv %*% (t(z) %*% mat %*% xk)
            xz[, k] <- xk + z %*% (dv_db[, k])
            dcla0b[, k] <- dcla0be %*% (xz[, k])
            dc <- Mi %*% diag(dcla0b[, k]) %*% done
            dmu0[, k] <- mu0 * (xz[, k])
            dw_db1[, k] <- dc * mu0 + mu * xz[, k]
            temp4 <- (2 * cla0/di) * (dcla0b[, k])
            dw_db2 <- ((diag(dmu0[, k])) %*% Mi %*% As %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% diag(temp4[, 1]) %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% As %*% t(Mi) %*% (diag(dmu0[, 
                k])))
            dw_db <- diag(dw_db1[, k]) - dw_db2
            ad1[k, 1] <- sum(diag(dinv %*% (t(z) %*% dw_db %*% 
                z)))
        }
        dft1 <- dft1 - 0.5 * ad1
    }
    dft2 <- crossprod(z, del - Wi %*% clam0) - (IArho %*% iD %*% v_h0)
    dft <- rbind(dft1, dft2)
    Bi <- diag(clam0[, 1])
    temp4 <- cla0^2/di
    Adi <- diag(temp4[, 1])
    As <- Adi
    mat <- (diag(Wi) * Bi) - (diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))
    U <- IArho%*%iD
#####################
    xmat<-crossprod(x,mat)
    zmat<-crossprod(z,mat)
    xmatx<-xmat%*%x
    xmatz<-xmat%*%z
    zmatx<-t(xmatz)
    zmatz<-zmat%*%z
    zmatz1<-zmat%*%z+U
    H <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz1))
    H0 <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz))
#####################
    Hinv <- solve(H)
    be_h0 <- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv %*% dft)
    beta_h[1:p, 1] <- be_h[1:p, 1]
    se_beta_h <- matrix(0, p, 1)
    for (i in 1:p) se_beta_h[i, 1] <- sqrt(Hinv[i, i])
    index2 <- qcum[nrand + 1]
    index3 <- p + 1
    index4 <- p + qcum[nrand + 1]
    v_h[1:index2, 1] <- be_h[index3:index4, 1]
    nonneg_adj=0
    if (varnonneg == TRUE) nonneg_adj=2
    for (i in 1:nrand) {
        if (dord == 0) {
            index1 <- p + qcum[i] + 1
            index2 <- p + qcum[i + 1]
            gamma <- sum(diag(Hinv[index1:index2, index1:index2]))/alpha_h[i]
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            if (varfixed == FALSE) 
                alpha_h[i] <- sum(v_h[index1:index2, 1]^2)/(q[i] - 
                  gamma-nonneg_adj)
        }
        if (dord == 1 | dord == 2) {
            H22 <- solve(zmatz1)
            ial1 <- 1/alpha_h[i]
            iA <- iD
            C <- matrix(0, qcum[nrand + 1], qcum[nrand + 1])
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            for (j in index1:index2) C[j, j] <- 1
            iB1 <- iA %*% C %*% iA
            c_vh1 <- iB1 %*% IArho %*% v_h
            dv1 <- H22 %*% c_vh1
            dexpeta1 <- expeta * (z %*% dv1)
            dcla01 <- -(di/(crossprod(Mi,expeta)^2)) * (crossprod(Mi,dexpeta1))
            dWi1 <- diag(dexpeta1[, 1])
            dAi1 <- diag(dcla01[, 1])
            temp4 <- Mi %*% dAi1 %*% done
            dBi1 <- diag(temp4[, 1])
            dvec1 <- 2 * (cla0 * dcla01)
            temp4 <- dvec1/di
            dAs1 <- diag(temp4[, 1])
            dmat1 <- (diag(dWi1) * Bi) + (diag(Wi) * dBi1) - (diag(dWi1)*Mi%*% (diag(As)*(t(diag(Wi)*Mi))))- (diag(Wi)*Mi %*% (diag(dAs1)*(t(diag(Wi)*Mi))))- (diag(Wi) * Mi %*% (diag(As)*(t(diag(dWi1)*Mi))))
            dia1 <- -iB1 %*% IArho
#####################
    xmat1<-crossprod(x,dmat1)
    zmat1<-crossprod(z,dmat1)
    xmat1x<-xmat1%*%x
    xmat1z<-xmat1%*%z
    zmat1x<-t(xmat1z)
    zmat1z<-zmat1%*%z
    zmat1z1<-zmat1%*%z+dia1
#####################
            Hd1 <- rbind(cbind(xmat1x, xmat1z), cbind(zmat1x, zmat1z1))
            gamma1 <- -alpha_h[i] * sum(Hinv * Hd1)
            if (dord == 2) {
                expeta <- exp(x %*% beta_h + z %*% v_h)
                ial1 <- 1/alpha_h[i]
                muu <- exp(x %*% beta_h) * clam0
                zmuu <- t(z) %*% muu
                u_h <- exp(v_h)
                ude1 <- u_h * dv1
                aa1 <- (zmuu * u_h) + ial1
                bb1 <- (zmuu * ude1) - (ial1^2)
                term11 <- ((aa1 * zmuu * ude1) - (2 * zmuu * 
                  u_h * bb1))/(aa1^3)
                term21 <- ((2 * aa1 * zmuu * zmuu * u_h * ude1) - 
                  (3 * ((zmuu * u_h)^2) * bb1))/(aa1^4)
                term1 <- (3 * term11) - (5 * term21)
                SS1 <- diag(term1[, 1])
                gamma21 <- -(alpha_h[i]/12) * sum(diag(SS1))
            }
            if (dord == 1) {
                gamma21 <- 0
            }
           k21 <- q[i] - gamma1 - gamma21 - t(v_h) %*% IArho %*% v_h /(alpha_h[i])
            if (varfixed == FALSE) 
                alpha_h[i] <- t(v_h) %*% IArho %*% v_h /(q[i] - 
                  gamma1 - gamma21-nonneg_adj)
        }
    }
    cumqrep <- cumsum(qrep)
    qq <- length(qrep)
    L1 <- 0
    L2 <- 0
    L3 <- 0
    n <- sum(qrep)
    for (i in 1:qq) {
        if (i==1) temp1 <- 1
        else temp1 <- cumqrep[i-1]+1
	temp2 <- cumqrep[i]
	Lambdai <- v_h[temp1:temp2,1] %*% t(v_h[temp1:temp2,1])
        temp3 <- nrow(beta_h)+temp1
        temp4 <- nrow(beta_h)+temp2
        Ti <- Hinv[temp3:temp4,temp3:temp4]
        Ii <- diag(1,qrep[i],qrep[i])
        Ji <- matrix(0,qrep[i],qrep[i])
        Ki <- matrix(0,qrep[i],qrep[i])
        for (k in 1:qrep[i]) {
            if (k<qrep[i]) Ji[k,k+1] <- 1
            if (k<qrep[i]) Ji[k+1,k] <- 1
            if (k==1) Ki[k,k] <- 1
            if (k==qrep[i]) Ki[k,k] <- 1
        }
        if (qrep[i]==1) Ki[1,1] <- 2
        L1 <- L1 + tr(Ii %*% (Ti + Lambdai))
        L2 <- L2 + 0.5*tr(Ji %*% (Ti + Lambdai))
        L3 <- L3 + tr(Ki %*% (Ti + Lambdai))
    }
    c1 <- (n-qq) * (L1-L3)
    c2 <- (2*qq - n) * L2
    c3 <- n*L3 - (n+qq)*L1
    c4 <- n*L2
    df_rho <- c1*rho^3 + c2*rho^2 + c3*rho + c4
    d2f_rho <- 3*c1*rho^2 + 2*c2*rho + c3
    rho <- rho - df_rho/d2f_rho
    se_rho <- sqrt(abs(-1/d2f_rho))
    res <- list(x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, 
        beta_h, v_h, alpha_h0, alpha_h, dft, Hinv, clam0, H, 
        mat, se_beta_h, U, H0, rho_0, rho, se_rho)
    return(res)
}
