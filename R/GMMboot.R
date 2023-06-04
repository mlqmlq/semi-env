dist <- function(A, B) {
  return(det(t(B) %*% A %*% t(A) %*% B))
}

#' Compute the mean of squared vector correlation among enveloped central space of bootstrap samples estimated by GMM method
#' @param X_data  predictor matrix
#' @param Y_data  Response matrix
#' @param u  Dimension of envelope selected
#' @param B  Number of bootstraps
#' @export
GMMboot<-function(X_data, Y_data, u, B){
  Y_data=as.matrix(Y_data)
  n = nrow(X_data)
  p = ncol(X_data)
  Gamma_u = nGMMenv(X_data, Y_data, u)
  Gamma_boot=list()
  q_2 = 0
  for (i in 1:B){
    resample = sample(1:n,n,replace=TRUE)
    X_data_r = X_data[resample,,drop=FALSE]
    Y_data_r = Y_data[resample,,drop=FALSE]
    Gamma_boot[[i]] = nGMMenv(X_data_r,Y_data_r,u)
    q_2 = q_2 + dist(Gamma_u$Gamma1, Gamma_boot[[i]]$Gamma1)
  }
  return(q_2/B)
}


#' Estimate the dimension of Enveloped central space via GMM estimator
#'
#' @param X  predictor matrix
#' @param Y  Response matrix
#' @param B  Number of bootstraps
#' @export
GMMdim = function(X, Y, B) {
  p=ncol(X)
  result=sapply(1:(p-1),function(t){
    GMMboot(X, Y, t, B=B)
  })
  list(dimsl=order(result)[p-1],VecCor=result)
}





