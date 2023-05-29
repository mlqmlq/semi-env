
# semienv


## Description

The goal of `semienv` is to derive
the entire class of regular and asymptotically linear (RAL) estimators as well as the
locally and globally semiparametrically efficient estimators for the enveloped central space. The package implements:

  * `nGMMenv3()`: GMM estimator; 
  * `localeff()`: semiparametrically local efficient estimator;
  * `globaleff()`: semiparametrically global efficient estimator
of enveloped central space. 

Additionally, this package also implements

  * `GMMdim()`: estimator of the dimension of Enveloped central space via the GMM estimator;
  * `effSDR()`: semiparametrically efficient estimator of dimension reduction subspace in Ma and Zhu ([2013](#semi-eff)).



## Installation
You can install the package with:

``` r
install.packages('devtools')
devtools::install_github("mlqmlq/semi-env")
```

## Usage
To demonstrate the usage of `semienv`  package, we provide several examples to show the procedures of estimating the enveloped central space via it. 

Given the simulated or real data, the first step of the estimation is determining the dimension of the enveloped central space by using `GMMdim()` function.

``` r
set.seed(500)
X=mvrnorm(500, mu = rep(0, 3), Sigma =diag(c(1,1,0.01)) , tol=1e-6)
Y=rbinom(500,1,1/(1+exp(-0.5*X[,1]-1*X[,2])))*(X[,1]+2*X[,2])
#Dimension Selection
u=GMMdim(X=X, ##Predictor
         Y=Y, ##Response          
         B=20 ##Number of bootstraps
         )$dimsl 
```
`GMMdim()` is a list that contains selected dimension `GMMdim()$dimsl` as well as sample mean of squared
vector correlation at each dimension `GMMdim()$VecCor`.

Then we can directly estimate enveloped central space via `nGMMenv3()`, `localeff()`, and `globaleff()` respectively with:
``` r
set.seed(500)
#GMM estimator
GMM_estimator=nGMMenv3(X=X, ##Predictor
                       Y=Y, ##Response
                       u=u ##selected dimension
                      )$Gamma1
#local efficient estimator
local_eff_estimator=localeff(X=X, ##Predictor
                       Y=Y, ##Response
                       beta=rnorm(3), #Model specification factor
                       u=u ##selected dimension
                            )$Gamma1  
#global efficient estimator
global_eff_estimator=globaleff(X=X, ##Predictor
                       Y=Y, ##Response                       
                       u=u ##selected dimension
                            )$Gamma1                                                      
```
All `nGMMenv3()`, `localeff()`, and `globaleff()` are lists which contain the  variational
independent parameter , `gamma`, and the semi-orthogonal matrix, `Gamma1`  of enveloped central space.

For more details regarding examples and arguments, please review the help page:

``` r
?semienv::nGMMenv3
?semienv::localeff
?semienv::globaleff
```

  

## Maintainer information

XXX XXXX ([XX@XX.edu](mailto:XX@XX.edu))


## References

<div id="refs" class="references">

<div id="semi-eff">

Ma, Y., and  Zhu, L. (2013). "Efficient estimation in sufficient dimension reduction." *The Annals of Statistics*, 41, 250-268.
<https://doi.org/10.1214/12-AOS1072>.

</div>
