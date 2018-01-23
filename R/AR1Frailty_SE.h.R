AR1Frailty_SE.h <-
function (res1, nrand, q, qcum, dord = 1, varfixed = FALSE, IArho) 
{
    x <- res1[1][[1]]
    z <- res1[2][[1]]
    y <- res1[3][[1]]
    del <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- res1[9][[1]]
    v_h <- res1[10][[1]]
    beta_h1 <- res1[11][[1]]
    v_h1 <- res1[12][[1]]
    alpha_h <- res1[13][[1]]
    alpha_h1 <- res1[14][[1]]
    dft <- res1[15][[1]]
    Hinv <- res1[16][[1]]
    clam0 <- res1[17][[1]]
    H <- res1[18][[1]]
    mat <- res1[19][[1]]
    U <- res1[21][[1]]
    H0 <- res1[22][[1]]
    n <- nrow(x)
    p <- ncol(x)
    u_h1 <- exp(v_h1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h1[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    iA <- iD
    Bi <- diag(clam0[, 1])
    muh <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(muh)
    cla0 <- di/(t(Mi) %*% expeta)
    temp4 <- cla0^2/di
    As <- diag(temp4[, 1])
    Wi <- diag(expeta[, 1])
    done <- matrix(1, idx2, 1)
    H22 <- solve(t(z) %*% mat %*% z + U)
    Hessian <- matrix(0, nrand, nrand)
    if (varfixed == FALSE) {
        se_lam <- rep(0, nrand)
        for (i in 1:nrand) se_lam[i] <- 0.000
    }
    if (varfixed == TRUE) {
        se_lam <- rep(0, nrand)
        for (i in 1:nrand) se_lam[i] <- 0.000
    }
    eta <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(eta)
    one <- matrix(1, n, 1)
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    pi <- 3.14159265359
    term0 <- t(Mi) %*% expeta
    hlike1 <- (t(one) %*% (del * eta)) - (t(done) %*% (di * log(term0)))
    hlike2 <- 0
    hlike3 <- 0
    for (i in 1:nrand) {
        if (alpha_h[i] > 1e-05) 
            hlike2 <- hlike2 - (q[i]/2) * log(2 * pi) - ((q[i]/2) * 
                log(alpha_h1[i])) + 0.5 * log(abs(det(IArho)))
        index1 <- qcum[i] + 1
        index2 <- qcum[i + 1]
        vv_h1 <- matrix(0, q[i], 1)
        vv_h1[1:q[i], 1] <- v_h1[index1:index2, 1]
        if (alpha_h[i] > 1e-05) 
            hlike3 <- hlike3 - (t(vv_h1) %*% IArho %*% vv_h1)/(2 * alpha_h1[i])
    }
    hliken <- hlike1 + hlike2 + hlike3
    cc1 <- svd(2 * pi * Hinv)
    for (i in 1:length(cc1$d)) if (cc1$d[i] < 1e-05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    adj1 <- 0.5 * logdet1
    hpn1 <- hliken + adj1
    muu <- exp(x %*% beta_h1) * clam0
    zmu <- t(z) %*% muu
    u_h1 <- exp(v_h1)
    second <- 0
    for (i in 1:nrand) {
        ialph1 <- 1/alpha_h1[i]
        a21 <- (zmu * u_h1) + ialph1
        b31 <- zmu * u_h1
        S11 <- 3 * (b31/(a21^2))
        S21 <- 5 * ((b31^2)/(a21^3))
        temp4 <- S11 - S21
        S31 <- diag(temp4[, 1])
        second <- second - sum(diag(S31))/24
    }
    H22 <- t(z) %*% mat %*% z + iD
    cc1 <- svd(H22/(2 * pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i] > 1e+05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    hpn2 <- hliken - 0.5 * logdet1
    hpn3 <- hpn1 + second
    df1 <- sum(diag(Hinv %*% H0))
    res <- list(se_lam, hliken, hpn1, hpn2, hpn3, hlike1, df1)
    return(res)
}
