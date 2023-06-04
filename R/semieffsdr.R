#' Estimate the dimension reduction direction via semiparametric efficient estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param ndr  Dimension of directions selected
#' @param beta_init Initial value of beta
#' @export
#' @references
#' Ma, Y., & Zhu, L. (2012). A semiparametric approach to dimension reduction. Journal of the American Statistical Association, 107(497), 168-179.
#'
#' Ma, Y., & Zhu, L. (2013). Efficient estimation in sufficient dimension reduction. Annals of statistics, 41(1), 250.
#' @examples
#' set.seed(123)
#' n=500
#' X=cbind(rnorm(n,sd=10),rnorm(n,sd=10),rnorm(n,sd=1),rnorm(n,sd=1))
#' Y=I(X[,1]+X[,2]>0)+rnorm(n,sd=0.5)
#' beta_init= orthoDr::initB(X, as.matrix(Y), ndr=1, bw = 1, method = "sir", ncore=0)
#' effsdr(X,Y,1,beta_init)
effsdr = function(X, Y, ndr, beta_init){
  p=ncol(X)
  n=nrow(X)
  dat1 = cbind(X %*% beta_init, Y)
  dat2 = cbind(X %*% beta_init)
  H10 = Hpi(dat1)
  H11 = Hpi(dat1, deriv.order = 1)
  bw = npregbw(xdat = X %*% beta_init, ydat = as.vector(Y))$bw

  if(ndr==1){
    # H20 = Hpi(data2)
    #H21 = Hpi(data2, deriv.order = 1)
    H20=H21=bw

  }

  if(ndr>1){
    H20 = Hpi(dat2)
    H21 = Hpi(dat2, deriv.order = 1)
  }

  getgamma<-function(Gamma){

    Gammarepr=Gamma%*%solve(Gamma[1:ncol(Gamma),])
    gammaM=Gammarepr[-(1:ncol(Gamma)),]
    gamma=as.vector(gammaM)
    return(gamma)
  }
  beta_init_f=(getgamma(beta_init))
  beta_init_fs=2/pi*atan(beta_init_f)

  semisdr2 = function(gamma) {

    gamma3=tan(pi/2*gamma)

    Beta1 = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
    #Beta1=beta
    # Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (ndr+1):p]
    data1 = cbind(X %*% Beta1, Y)
    data2 = X %*% Beta1

    den1 = kde(data1, H10, eval.points = data1)$estimate
    den2 = kde(data2, H10[1:ndr,1:ndr], eval.points = data2)$estimate

    # denbx=kde(data2, H20, eval.points = data2,w=data2*length(data2)/sum(data2))$estimate
    conden=den1/den2

    dconden = t(t(kdde(data1, H10, deriv.order = 1, eval.points = data1)$estimate[,1:ndr]) * (den2)-t(kdde(data2, H10[1:ndr,1:ndr], deriv.order = 1, eval.points = data2)$estimate)*(den1))* (1/den2**2)
    dlogcd=dconden*(1/conden)


    QB=as.matrix(diag(p))- as.matrix(Beta1 %*% ginv(t(Beta1) %*% Beta1) %*% t(Beta1))
    QX = X %*%QB
    for (i in 1:p){
      denbx = kde(data2, H20, eval.points = data2,w=X[,i]*length(X[,i])/sum(X[,i]))$estimate
      QX[,i]=X[,i]-denbx*sum(X[,i])/(length(X[,i])*den2)
    }

    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogcd))] * dlogcd[, rep(seq(ncol(dlogcd)),ncol(QX))]

    return(-sum(log(den1/den2)))
  }
  gamma_prs= optim(par =beta_init_fs, fn = semisdr2 , method = "L-BFGS-B", lower = rep(-0.9999999, ndr*(p-ndr)), upper =  rep(0.999999, ndr*(p-ndr)))$par

  gamma_pr= tan(pi/2*gamma_prs)
  Beta_pr = qr.Q(qr(rbind(diag(ndr), matrix(gamma_pr, ncol = ndr))))
  Beta0_pr = qr.Q(qr(Beta_pr), complete = T)[, (ndr+1):p]
  Xg = cbind(X %*% Beta_pr,X %*% Beta0_pr)

  dat1 = cbind(Xg[,1:ndr], Y)
  dat2 = cbind(Xg[,1:ndr])
  H10 = Hpi(dat1)
  H11 = Hpi(dat1, deriv.order = 1)
  bw = npregbw(xdat = Xg %*% Beta_pr, ydat = as.vector(Y))$bw
  if(ndr==1){
    # H20 = Hpi(data2)
    #H21 = Hpi(data2, deriv.order = 1)
    H20=H21=bw

  }

  if(ndr>1){
    H20 = Hpi(dat2)
    H21 = Hpi(dat2, deriv.order = 1)
  }




  semisdr3 = function(gamma) {

    gamma3=tan(pi/2*gamma)

    Beta1 = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
    #Beta1=beta
    # Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (ndr+1):p]
    data1 = cbind(Xg %*% Beta1, Y)
    data2 = Xg %*% Beta1

    den1 = kde(data1, H10, eval.points = data1)$estimate
    den2 = kde(data2, H10[1:ndr,1:ndr], eval.points = data2)$estimate
    den1d=matrix(rep(0,(ndr)*n),nrow =n)
    den2d=matrix(rep(0,ndr*n),nrow=n)
    # denbx=kde(data2, H20, eval.points = data2,w=data2*length(data2)/sum(data2))$estimate
    conden=den1/den2

    dconden = t(t(kdde(data1, H10, deriv.order = 1, eval.points = data1)$estimate[,1:ndr]) * (den2)-t(kdde(data2, H10[1:ndr,1:ndr], deriv.order = 1, eval.points = data2)$estimate)*(den1))* (1/den2**2)
    dlogcd=dconden*(1/conden)


    dlogcd2=den1d/den1-den2d/den2


    QB=as.matrix(diag(p))- as.matrix(Beta1 %*% ginv(t(Beta1) %*% Beta1) %*% t(Beta1))
    QX = Xg %*%QB
    for (i in 1:p){
      denbx = kde(data2, H20, eval.points = data2,w=Xg[,i]*length(Xg[,i])/sum(Xg[,i]))$estimate
      QX[,i]=Xg[,i]-denbx*sum(Xg[,i])/(length(Xg[,i])*den2)
    }

    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogcd))] * dlogcd[, rep(seq(ncol(dlogcd)),ncol(QX))]

    crf=solve((-gamma3%*%t(gamma3)+diag(diag(gamma3%*%t(gamma3)))+diag((1+sum(gamma3**2)-gamma3**2)))/(1+sum(gamma3**2))**(3/2))
    return( (sum(abs(apply(vec1%*%dvec(gamma3, ndr)%*%crf, 2, mean)))))
  }
  gamma = optim(par =rep(0, ndr*(p-ndr)), fn = semisdr3 , method = "L-BFGS-B", lower = rep(-0.999, ndr*(p-ndr)), upper =  rep(0.999, ndr*(p-ndr)))$par
  gamma3=tan(pi/2*gamma)
  beta = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
  list(beta = (cbind(Beta_pr,Beta0_pr))%*%beta)
}










effsdr2 = function(X, Y, ndr, beta_init){
  p=ncol(X)
  n=nrow(X)
  dat1 = cbind(X %*% beta_init, Y)
  dat2 = cbind(X %*% beta_init)
  H10 = Hpi(dat1)
  H11 = Hpi(dat1, deriv.order = 1)
  bw = npregbw(xdat = X %*% beta_init, ydat = as.vector(Y))$bw

  if(ndr==1){
    # H20 = Hpi(data2)
    #H21 = Hpi(data2, deriv.order = 1)
    H20=H21=bw

  }

  if(ndr>1){
    H20 = Hpi(dat2)
    H21 = Hpi(dat2, deriv.order = 1)
  }

  getgamma<-function(Gamma){

    Gammarepr=Gamma%*%solve(Gamma[1:ncol(Gamma),])
    gammaM=Gammarepr[-(1:ncol(Gamma)),]
    gamma=as.vector(gammaM)
    return(gamma)
  }
  beta_init_f=(getgamma(beta_init))
  beta_init_fs=2/pi*atan(beta_init_f)

  semisdr2 = function(gamma) {

    gamma3=tan(pi/2*gamma)

    Beta1 = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
    #Beta1=beta
    # Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (ndr+1):p]
    data1 = cbind(X %*% Beta1, Y)
    data2 = X %*% Beta1

    den1 = kde(data1, H10, eval.points = data1)$estimate
    den2 = kde(data2, H10[1:ndr,1:ndr], eval.points = data2)$estimate

    # denbx=kde(data2, H20, eval.points = data2,w=data2*length(data2)/sum(data2))$estimate
    conden=den1/den2

    dconden = t(t(kdde(data1, H10, deriv.order = 1, eval.points = data1)$estimate[,1:ndr]) * (den2)-t(kdde(data2, H10[1:ndr,1:ndr], deriv.order = 1, eval.points = data2)$estimate)*(den1))* (1/den2**2)
    dlogcd=dconden*(1/conden)


    QB=as.matrix(diag(p))- as.matrix(Beta1 %*% ginv(t(Beta1) %*% Beta1) %*% t(Beta1))
    QX = X %*%QB
    for (i in 1:p){
      denbx = kde(data2, H20, eval.points = data2,w=X[,i]*length(X[,i])/sum(X[,i]))$estimate
      QX[,i]=X[,i]-denbx*sum(X[,i])/(length(X[,i])*den2)
    }
    #  return(log(abs(sum(dlogcd))))
    # QX = X %*%QB
    # denbx = kde(data2, H20, eval.points = data2,w=data2*n/sum(data2))$estimate
    #QX = X %*%Q(Beta1)
    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogcd))] * dlogcd[, rep(seq(ncol(dlogcd)),ncol(QX))]
    # vec1 = X[, rep(seq(ncol(X)), each = ncol(dlogcd))] * dlogcd2[, rep(seq(ncol(dlogcd)),ncol(X))]

    #dlogden2 = t(t(kdde(data2, H21, deriv.order = 1, eval.points = data2)$estimate) * (1/den2))
    # PX = X %*% P(Gamma1)
    #  vec2 = PX[, rep(seq(ncol(PX)), each = ncol(dlogden2))] * dlogden2[, rep(seq(ncol(dlogden2)),ncol(PX))]
    #return( sum((cov(QX,dlogcd))))
    return(-sum(log(den1/den2)))
    #return( (sum(abs(apply(vec1%*%dvec(gamma3, ndr)*(1/cos(pi/2*gamma)**2), 2, mean)))))
  }
  ga= optim(par =beta_init_fs, fn = semisdr2 , method = "L-BFGS-B", lower = rep(-0.9999999, ndr*(p-ndr)), upper =  rep(0.999999, ndr*(p-ndr)))$par


  semisdr3 = function(gamma) {

    gamma3=tan(pi/2*gamma)

    Beta1 = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
    #Beta1=beta
    # Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (ndr+1):p]
    data1 = cbind(X %*% Beta1, Y)
    data2 = X %*% Beta1

    den1 = kde(data1, H10, eval.points = data1)$estimate
    den2 = kde(data2, H10[1:ndr,1:ndr], eval.points = data2)$estimate
    den1d=matrix(rep(0,(ndr)*n),nrow =n)
    den2d=matrix(rep(0,ndr*n),nrow=n)
    # denbx=kde(data2, H20, eval.points = data2,w=data2*length(data2)/sum(data2))$estimate
    conden=den1/den2

    dconden = t(t(kdde(data1, H10, deriv.order = 1, eval.points = data1)$estimate[,1:ndr]) * (den2)-t(kdde(data2, H10[1:ndr,1:ndr], deriv.order = 1, eval.points = data2)$estimate)*(den1))* (1/den2**2)
    dlogcd=dconden*(1/conden)


    dlogcd2=den1d/den1-den2d/den2


    QB=as.matrix(diag(p))- as.matrix(Beta1 %*% ginv(t(Beta1) %*% Beta1) %*% t(Beta1))
    QX = X %*%QB
    for (i in 1:p){
      denbx = kde(data2, H20, eval.points = data2,w=X[,i]*length(X[,i])/sum(X[,i]))$estimate
      QX[,i]=X[,i]-denbx*sum(X[,i])/(length(X[,i])*den2)
    }
    #  return(log(abs(sum(dlogcd))))
    # QX = X %*%QB
    # denbx = kde(data2, H20, eval.points = data2,w=data2*n/sum(data2))$estimate
    #QX = X %*%Q(Beta1)
    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogcd))] * dlogcd[, rep(seq(ncol(dlogcd)),ncol(QX))]
    # vec1 = X[, rep(seq(ncol(X)), each = ncol(dlogcd))] * dlogcd2[, rep(seq(ncol(dlogcd)),ncol(X))]

    #dlogden2 = t(t(kdde(data2, H21, deriv.order = 1, eval.points = data2)$estimate) * (1/den2))
    # PX = X %*% P(Gamma1)
    #  vec2 = PX[, rep(seq(ncol(PX)), each = ncol(dlogden2))] * dlogden2[, rep(seq(ncol(dlogden2)),ncol(PX))]
    #return( sum((cov(QX,dlogcd))))
    #return(-sum(log(den1/den2)))
    return( (sum(abs(apply(vec1%*%dvec(gamma3, ndr)*(1/cos(pi/2*gamma)**2), 2, mean)))))
  }
  gamma = optim(par =ga, fn = semisdr3 , method = "L-BFGS-B", lower = rep(-0.9999999, ndr*(p-ndr)), upper =  rep(0.999999, ndr*(p-ndr)))$par
  gamma3=tan(pi/2*gamma)
  beta = qr.Q(qr(rbind(diag(ndr), matrix(gamma3, ncol = ndr))))
  list(beta = beta)
}

