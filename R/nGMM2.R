#' Estimate enveloped central space via GMM estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param u  Dimension of envelope selected
#' @export
#' @examples
#' #generate  data
#' set.seed(1)
#' Gamma = diag(3)[,1:2]
#' Gamma0 = diag(3)[,3]
#' Omega = matrix(1, nrow = 2, ncol = 2) + diag(2, 2)
#' Omega0 = diag(0.1, 1)
#' SigmaX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
#' X=mvrnorm(100, mu = rep(0, 3), Sigma = SigmaX, tol=1e-6)
#' Y=X[,1]+rnorm(100)
#' Get envelope space
#' result=nGMMenv(X,Y,2)
#' Gammaest=result$Gamma1
#' Gammaest
nGMMenv = function(X, Y, u) {
  p=ncol(X)
  gmm_obj = function(gamma) {
    Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u))))
    Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (u+1):p]
    sum(abs(c(c(cor(Y, X %*% Gamma0)), c(cor(X %*% Gamma1, X %*% Gamma0)))))
  }
  Gamma_cdd = vector('list', p)
  reodd = numeric()
  Gamma_cdd[[1]] = optim(par = rep(0, u*(p-u)), fn = gmm_obj, method = "BFGS")$par
  reodd[1] = gmm_obj(Gamma_cdd[[1]])
  for(i in 2:(2*p+1)) {
    Gamma_cdd[[i]] = optim(par = runif(u*(p-u), -0.25, 0.25), fn = gmm_obj, method = "BFGS")$par
    reodd[i] = gmm_obj(Gamma_cdd[[i]])
  }
  gamma = Gamma_cdd[[order(reodd)[1]]]
  list(gamma = gamma, Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u)))))
}

nGMMenv2 = function(X, Y, u) {
  p=ncol(X)
  gmm_obj3 = function(gamma) {
    gamma3=tan(pi/2*gamma)
    Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma3, ncol = u))))
    Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (u+1):p]
    sum(abs(c(c(cor(Y, X %*% Gamma0)), c(cor(X %*% Gamma1, X %*% Gamma0)))))
  }
  Gamma_cdd = vector('list', p)
  reodd = numeric()
  Gamma_cdd[[1]] = optim(par = rep(0, u*(p-u)), fn = gmm_obj3, method = "L-BFGS-B", lower = rep(-1, u*(p-u)), upper = rep(1, u*(p-u)))$par
  reodd[1] = gmm_obj3(Gamma_cdd[[1]])
  for(i in 2:(p**2+1)) {
    Gamma_cdd[[i]] = optim(par = runif(u*(p-u), -1, 1), fn = gmm_obj3, method = "L-BFGS-B", lower = rep(-1, u*(p-u)), upper = rep(1, u*(p-u)))$par
    reodd[i] = gmm_obj3(Gamma_cdd[[i]])
  }
  gamma = tan(pi/2*Gamma_cdd[[order(reodd)[1]]])
  list(gamma = gamma, Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u)))))
}




#' Estimate enveloped central space via GMM estimator, tangent scaling parameter
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param u  Dimension of envelope selected
#' @export
#' @examples
#' #generate  data
#' set.seed(3)
#' Gamma = diag(5)[,c(1,4)]
#' Gamma0 = diag(5)[,c(2,3,5)]
#' Omega = matrix(1, nrow = 2, ncol = 2) + diag(2, 2)
#' Omega0 = diag(c(0.2,0.2,2))
#' SigmaX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
#' X=mvrnorm(500, mu = rep(0, 5), Sigma = SigmaX, tol=1e-6)
#' Y=exp(X[,1])+8*I(X[,3]>0)+rnorm(500)
#' #Get Envelope Space
#' result=nGMMenv3(X,Y,3)
#' Gammaest=result$Gamma1
#' Gammaest
nGMMenv3 = function(X, Y, u) {
  p=ncol(X)
  gmm_obj3 = function(gamma) {
    gamma3=tan(pi/2*gamma)
    Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma3, ncol = u))))
    Gamma0 = qr.Q(qr(Gamma1), complete = T)[, (u+1):p]
    sum(abs(c(c(cor(Y, X %*% Gamma0)), c(cor(X %*% Gamma1, X %*% Gamma0)))))
  }
  Gamma_cdd = vector('list', p)
  reodd = numeric()
  Gamma_cdd[[1]] = optim(par = rep(0, u*(p-u)), fn = gmm_obj3, method = "L-BFGS-B", lower = rep(-1, u*(p-u)), upper = rep(1, u*(p-u)))$par
  reodd[1] = gmm_obj3(Gamma_cdd[[1]])
  for(i in 2:(p**2+1)) {
    Gamma_cdd[[i]] = optim(par = runif(u*(p-u), -1, 1), fn = gmm_obj3, method = "L-BFGS-B", lower = rep(-1, u*(p-u)), upper = rep(1, u*(p-u)))$par
    reodd[i] = gmm_obj3(Gamma_cdd[[i]])
  }
  gamma = tan(pi/2*Gamma_cdd[[order(reodd)[1]]])
  list(gamma = gamma, Gamma1 = qr.Q(qr(rbind(diag(u), matrix(gamma, ncol = u)))))
}






