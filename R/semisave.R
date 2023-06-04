#' Estimate the dimension reduction direction via semi-SAVE estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param u  Dimension of directions selected
#' @param beta_init Initial value of beta
#' @export
semisave = function(X, Y, ndr, beta_init){
  dat1 = cbind(X %*% beta_init, Y)
  dat2 = cbind(X %*% beta_init)
  H10 = Hpi(dat1)
  H11 = Hpi(dat1, deriv.order = 1)
  bw = npregbw(xdat = X %*% beta_init, ydat = as.vector(Y))$bw
  
  if(ndr==1){
    # H20 = Hpi(data2)
    #H21 = Hpi(data2, deriv.order = 1)
    H20=H21=bw
    bwy=bw

  }
  
  if(ndr>1){
    H20 = Hpi(dat2)
    H21 = Hpi(dat2, deriv.order = 1)
    bwy=hpi(Y)
  }
  
  
  
 
  savef3 = function(beta) {
   # Beta1 = qr.Q(qr(rbind(diag(ndr), matrix(gamma, ncol = ndr))))
     Beta1=matrix(as.vector(beta),nrow=p)
    # Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (ndr+1):p]
    #data1 = cbind(X %*% Beta1, Y)
    data2 = X %*% Beta1
    #  den1 = kde(data1, H10, eval.points = data1)$estimate
    den2 = kde(data2, H20, eval.points = data2)$estimate
    deny = kde(Y, bw, eval.points = Y)$estimate
    conden=den1/den2
    #dconden = t(t(kdde(data1, H11, deriv.order = 1, eval.points = data1)$estimate[,1:ndr]) * (den2)-t(kdde(data2, H21, deriv.order = 1, eval.points = data2)$estimate)*(den1))* (1/den2**2)
    #dlogcd=dconden*(1/conden)
    #QB=as.matrix(diag(p))- as.matrix(Beta1 %*% ginv(t(Beta1) %*% Beta1) %*% t(Beta1))
    #QX = X %*%QB
    #QX = X %*%Q(Beta1)
    
    cdxxbx=matrix(rep(0,p**2*n),ncol=p**2)
    cdxxy=matrix(rep(0,p**2*n),ncol=p**2)
    cdxbx=matrix(rep(0,p*n),ncol=p)
    cdxy=matrix(rep(0,p*n),ncol=p)
    
    for (i in 1:p){
      cdxbx[,i] = (kde(data2, H20, eval.points = data2,w=X[,i]*length(X[,i])/sum(X[,i]))$estimate)*sum(X[,i])/(length(X[,i])*den2)
      cdxy[,i] = (kde(Y, bwy, eval.points = Y,w=X[,i]*length(Y)/sum(X[,i]))$estimate)*sum(X[,i])/(length(Y)*deny)
      for (j in 1:p){
        cdxxy[,(i-1)*p+j]= (kde(Y, bwy, eval.points = Y,w=X[,i]*X[,j]*length(Y)/sum(X[,i]*X[,j]))$estimate)*sum(X[,i]*X[,j])/(length(Y)*deny)
        cdxxbx[,(i-1)*p+j]= (kde(data2, H20, eval.points = data2,w=X[,i]*X[,j]*length(X[,i])/sum(X[,i]*X[,j]))$estimate)*sum(X[,i]*X[,j])/(length(X[,i])*den2)
        
      }
    }
    
    vec2=matrix(rep(0,p**2*n),ncol=p**2)
    for (k in 1:length(X[,1])){
      m1=as.matrix(diag(p))-matrix(cdxxy[k,],ncol=p)+cdxy[k,]%*%t(cdxy[k,])
      m2=(X[k,]-cdxy[k,])%*%t((X[k,]-cdxbx[k,]))-matrix(cdxxbx[k,],ncol=p)+cdxbx[k,]%*%t(cdxbx[k,])
      vec2[k,]=as.vector(m1%*%m2)
      
    }
    
    
    
    
    
    
    
    # vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogcd))] * dlogcd[, rep(seq(ncol(dlogcd)),ncol(QX))] 
    
    #dlogden2 = t(t(kdde(data2, H21, deriv.order = 1, eval.points = data2)$estimate) * (1/den2))
    # PX = X %*% P(Gamma1)
    #  vec2 = PX[, rep(seq(ncol(PX)), each = ncol(dlogden2))] * dlogden2[, rep(seq(ncol(dlogden2)),ncol(PX))] 
    return( sum(abs(apply(vec2, 2, mean))))
  }
  beta = optim(par = as.vector(beta_init), fn = savef3 , method = "L-BFGS-B", lower = as.vector(beta_init + rep(-0.05, ndr*(p))), upper = as.vector(beta_init+ rep(0.05, ndr*(p))))$par
  beta= matrix(as.vector(beta),nrow=p)
  list(beta = beta)
}
