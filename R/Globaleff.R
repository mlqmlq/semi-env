#' Compute the projection matrix onto the column space of the design matrix A
#' @param A the design matrix A
#' @export
P = function(A) { A %*% ginv(t(A) %*% A) %*% t(A) }
#' Compute the projection matrix onto the orthogonal complement of the column space of the design matrix A
#' @param A the design matrix A
#' @export
Q = function(A) { diag(nrow(A)) - A %*% ginv(t(A) %*% A) %*% t(A) }


dvec = function(t, u) {
  h = 5e-6
  y = c(qr.Q(qr(rbind(diag(u), matrix(t, ncol = u)))))
  out = matrix(0, nrow = length(y), ncol = length(t))
  for (i in 1:length(t)) {
    t1 = t
    t1[i] = t1[i] + h
    y1 = c(qr.Q(qr(rbind(diag(u), matrix(t1, ncol = u)))))
    t2 = t
    t2[i] = t2[i] - h
    y2 = c(qr.Q(qr(rbind(diag(u), matrix(t2, ncol = u)))))
    out[, i] = (y1 - y2) / (2*h)
  }
  return(out)
}

dvec0 = function(t, u) {
  h = 5e-6
  y = c(qr.Q(qr(rbind(diag(u), matrix(t, ncol = u))), complete = T)[, -(1:u)])
  out = matrix(0, nrow = length(y), ncol = length(t))
  for (i in 1:length(t)) {
    t1 = t
    t1[i] = t1[i] + h
    y1 = c(qr.Q(qr(rbind(diag(u), matrix(t1, ncol = u))), complete = T)[, -(1:u)])
    t2 = t
    t2[i] = t2[i] - h
    y2 = c(qr.Q(qr(rbind(diag(u), matrix(t2, ncol = u))), complete = T)[, -(1:u)])
    out[, i] = (y1 - y2) / (2*h)
  }
  return(out)
}


#' Estimate enveloped central space via global efficient estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param gamma_init  Initial value of gamma, defaultly obtaining from GMM method
#' @param u  Dimension of Envelope of selected
#' @export
#' @examples
#' set.seed(500)
#' X=mvrnorm(500, mu = rep(0, 3), Sigma =diag(c(1,1,0.01)) , tol=1e-6)
#' Y=rbinom(500,1,1/(1+exp(-0.5*X[,1]-1*X[,2])))*(X[,1]+2*X[,2])
#' result=globaleff(X,Y,1)
#' Gammaest=result$Gamma1
#' Gammaest

globaleff = function(X, Y,u,gamma_init="notgiven"){
  n=nrow(X)
  C=1
  p=ncol(X)

  if((gamma_init[1]=="notgiven")|(!is.vector(gamma_init))|(length(gamma_init)!=u*(p-u))){
    print("initial gamma is not given or the format is wrong,we defaultly use GMM method as initial value")
    gamma_init=nGMMenv2(X,Y,u)$gamma

  }

  Gamma_init = qr.Q(qr(rbind(diag(u), matrix(gamma_init, ncol = u))))
  Gamma0_init = qr.Q(qr(Gamma_init), complete = T)[, (u+1):p]

  data1 = cbind(X %*% Gamma_init, Y)
  data2 = cbind(X %*% Gamma0_init)
  H10 = Hpi(data1)
  H11 = Hpi(data1, deriv.order = 1)
  if(p-u>1){
    H20 = Hpi(data2)
    H21 = Hpi(data2, deriv.order = 1)
  }
  if(p-u==1){
    H20=hpi(data2, deriv.order = 1)
    H21=hpi(data2, deriv.order = 1)
  }

  semiglobal = function(gamma) {
    Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u))))
    Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (u+1):p]
    data1 = cbind(X %*% Gamma1, Y)
    data2 = X %*% Gamma0
    den1 = kde(data1, H10, eval.points = data1)$estimate
    dlogden1 = t(t(kdde(data1, H11, deriv.order = 1, eval.points = data1)$estimate[,1:u]) * (1/den1))
    QX = X %*% Q(Gamma1)
    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlogden1))] * dlogden1[, rep(seq(ncol(dlogden1)),ncol(QX))]
    den2 = kde(data2, H20, eval.points = data2)$estimate
    dlogden2 = t(t(kdde(data2, H21, deriv.order = 1, eval.points = data2)$estimate) * (1/den2))
    PX = X %*% P(Gamma1)
    vec2 = PX[, rep(seq(ncol(PX)), each = ncol(dlogden2))] * dlogden2[, rep(seq(ncol(dlogden2)),ncol(PX))]
  # -sum(log(den1)+log(den2))
     sum(abs(apply(vec1 %*% dvec(gamma, u) + vec2 %*% dvec0(gamma, u), 2, mean)))

  }
  gamma = optim(par = gamma_init, fn = semiglobal, method = "L-BFGS-B", lower = gamma_init-rep(C/sqrt(n), u*(p-u)), upper =gamma_init+rep(C/sqrt(n), u*(p-u)))$par
  list(gamma = gamma, Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u)))))
}



#' Estimate enveloped central space via local efficient estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param gamma_init  Initial value of gamma, defaultly obtaining from GMM method
#' @param beta Model specification factor
#' @param u  Dimension of envelope selected
#' @export
#' @examples
#'set.seed(300)
#'Gamma0 = diag(4)[,c(1,4)]
#'Gamma = diag(4)[,c(2,3)]
#'Omega = matrix(1, nrow = 2, ncol = 2) + diag(2, 2)
#'Omega0 = diag(c(2,2))
#'SigmaX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
#'X=mvrnorm(300, mu = rep(0, 4), Sigma = SigmaX, tol=1e-6)
#'Y=numeric()
#'for(i in 1:300)
#'Y[i]=sample(c(0,X[i,3]/2,X[i,3]),1)+rnorm(1)
#'result=localeff(X,Y,2,beta=c(1,1,1,1))
#result$Gamma1
localeff = function(X, Y, u, gamma_init ="notgiven", beta) {
  n=nrow(X)
  C=10
  p=ncol(X)


  if((gamma_init[1]=="notgiven")|(!is.vector(gamma_init))|(length(gamma_init)!=u*(p-u))){
    print("initial gamma is not given or the format is wrong,we defaultly use GMM method as initial value")
    gamma_init=nGMMenv2(X,Y,u)$gamma

  }

  Gamma_init = qr.Q(qr(rbind(diag(u), matrix(gamma_init, ncol = u))))

  bw = npregbw(xdat = X %*% Gamma_init, ydat = as.vector(Y))$bw

  semilocal = function(gamma) {
    Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u))))
    Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (u+1):p]
    eta = t(Gamma1) %*% beta
    SigmaX = var(X)
    Omega0 = t(Gamma0) %*% SigmaX %*% Gamma0
    sigma2 = var(Y - X %*% beta)[1]
    np1 = npreg(xdat = X %*% Gamma1, ydat = as.vector(Y), bws = bw)
    dlog1 = kronecker(t(eta), Y - predict(np1)) / sigma2
    QX = X %*% Q(Gamma1)
    vec1 = QX[, rep(seq(ncol(QX)), each = ncol(dlog1))] * dlog1[, rep(seq(ncol(dlog1)),ncol(QX))]
    dlog2 = kronecker(t(apply(X %*% Gamma0 %*% ginv(Omega0), 2, mean)), matrix(1, nrow = nrow(X))) - X %*% Gamma0 %*% ginv(Omega0)
    PX = X %*% P(Gamma1)
    vec2 = PX[, rep(seq(ncol(PX)), each = ncol(dlog2))] * dlog2[, rep(seq(ncol(dlog2)),ncol(PX))]
    sum(abs(apply(vec1 %*% dvec(gamma, u) + vec2 %*% dvec0(gamma, u), 2, mean)))
  }
  gamma = optim(par = gamma_init, fn = semilocal, method = "L-BFGS-B", lower = gamma_init-rep(C/sqrt(n), u*(p-u)), upper =gamma_init+rep(C/sqrt(n), u*(p-u)))$par
  list(gamma = gamma, Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u)))))
}










