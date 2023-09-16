#' @title Genetic correlation class
#' @description This class is only defined for the estimation of genetic correlation, genetic covariance,
#'  heritability, genetic coefficient.
#' @param input_parameter <Input Parameter Description>
#' @export <Do you want users to be able to use this function? If not, it is a
#' @keywords
#' @seealso
#' @import <Load entire package as dependency>
#' @importFrom <Load specific functions from packages as dependencies>
#' @return Object returned
#' @aliases
#' @examples
#'  n.z = sum(!is.na(Dataset$Z))
#'  n.y = sum(!is.na(Dataset$Y))
#'  X = as.matrix(Dataset[,-c(1,2)])
#'  Y.hat <- f(X %*% beta.hat[-1] + beta.hat[1] )
#'  Z.hat <-  g(X %*%gamma.hat[-1] + gamma.hat[1])
#'  eps.hat <- Dataset$Y - Y.hat;v.hat <- Dataset$Z - Z.hat;
#'  eps.hat[is.na(eps.hat)] = 0;v.hat[is.na(v.hat)] = 0;
#'  GePa= GeCr$new(N.Y = n.y, N.Z = n.z, 
#'   Y.hat = Y.hat, 
#'   Z.hat = Z.hat, 
#'   eps.hat = eps.hat, 
#'   v.hat = v.hat)

library(R6)

##############################################
#
# Calculating the esimated variance of heritability
#
##############################################

GeCr <- R6Class("GeCr", list(
  N.Y = NA,
  N.Z = NA,
  N = NA,
  Y.hat = NA,
  Z.hat = NA,
  eps.hat = NA,
  v.hat = NA,
  Delta.h2.beta = NA,
  Delta.h2.gamma = NA,
  Delta.gecv = NA,
  h2.beta.est = NA,
  h2.gamma.est = NA,
  Q.beta.est = NA,
  Q.gamma.est = NA,
  GeCv.est = NA,
  GeCr.est = NA,
  IV.cef.est = NA,
  h2.beta.est.se = NA,
  h2.gamma.est.se = NA,
  Q.beta.est.se = NA,
  Q.gamma.est.se = NA,
  GeCv.est.se = NA,
  GeCr.est.se = NA,
  IV.cef.est.se = NA,
  weight = 1,
  initialize = function(N.Y, N.Z, Y.hat, Z.hat, eps.hat, v.hat, weight = 1){
    
    stopifnot(is.numeric(N.Y), is.numeric(N.Z))
    stopifnot( length(Y.hat) == length(Z.hat), length(eps.hat) == length(v.hat) )
    self$Y.hat = Y.hat; self$Z.hat = Z.hat; self$eps.hat =  eps.hat; self$v.hat =  v.hat;
    self$N.Y = N.Y; self$N.Z = N.Z;self$N = length(Y.hat);
    self$weight = weight
    mu.y.hat = mean(self$weight * self$Y.hat) + sum( self$weight * self$eps.hat)/self$N.Y
    mu.z.hat = mean(self$weight * self$Z.hat) + sum( self$weight * self$v.hat)/self$N.Z
    self$Delta.h2.beta = self$weight *(self$Y.hat - mu.y.hat) * (self$Y.hat - mu.y.hat) /self$N + 
      2*self$weight * self$eps.hat * (self$Y.hat - mu.y.hat) /self$N.Y
    
    self$Delta.h2.gamma = self$weight*(self$Z.hat - mu.z.hat) * (self$Z.hat - mu.z.hat) /self$N + 
      2*self$weight * self$v.hat * (self$Z.hat - mu.z.hat)/self$N.Z
    
    self$Delta.gecv = self$weight *(self$Y.hat - mu.y.hat)  * (self$Z.hat - mu.z.hat) /self$N + 
      self$weight * self$eps.hat * (self$Z.hat - mu.z.hat) /self$N.Y + 
      self$weight * self$v.hat * (self$Y.hat - mu.y.hat)/self$N.Z

    self$h2.beta.est = sum(self$Delta.h2.beta)
    self$h2.beta.est.se = sqrt(  length(self$Delta.h2.beta)*var(self$Delta.h2.beta ) + (sum(self$Delta.h2.beta) -  self$h2.beta.est)^2/self$N )
    self$h2.gamma.est = sum(self$Delta.h2.gamma)
    self$h2.gamma.est.se = sqrt(  length(self$Delta.h2.gamma)*var(self$Delta.h2.gamma) + (sum(self$Delta.h2.gamma) -  self$h2.gamma.est)^2/self$N )
    self$GeCv.est = sum(self$Delta.gecv)
    self$GeCv.est.se =  sqrt(  length(self$Delta.gecv)*var(self$Delta.gecv) + (sum(self$Delta.gecv) -  self$GeCv.est)^2/self$N )
  },
  print = function(...){
    cat(" A class for genetic correlation:\n")
    cat("The heritability  of outomce 1 is: ",  self$h2.beta.est, "\n", "standard error is ", self$h2.beta.est.se,"\n", sep = "")
    cat("The heritability  of outcome 2 is: ",  self$h2.gamma.est, "\n", "standard error is ", self$h2.gamma.est.se,"\n", sep = "")
    cat("The Genetic Covariance between outcome 1 and outcome 2 is: ",  self$GeCv.est, "\n", "standard error is ", self$GeCv.est.se, "\n", sep = "")
    cat("The Genetic Correlation between outcome 1 and outcome 2 is: ",  self$GeCr.est, "\n", "standard error is ", self$GeCr.est.se, sep = "")
  },
  cov.est.hat = function(Delta){
    #' @description Estimator for the GeCov
    
    estimator = sum(Delta)
    estimator.var = length(Delta)*var(Delta ) + (sum(Delta) -  estimator)^2/sum(!is.na(Delta))
    
    # plugest = self$plug() - mu.y*mu.z
    return(c(estimator, estimator.var))
  },
  GeCr.est.hat = function(){
    #' @description Estimator for the GeCr
    
    self$GeCr.est = self$GeCv.est/sqrt(self$h2.beta.est *self$h2.gamma.est)
    
    
    
    
    #--------------------- Jacknife estimate of the variance ---------------#
    delta.mat = data.frame(Delta.h2.beta = self$Delta.h2.beta, 
                           Delta.h2.gamma = self$Delta.h2.gamma, 
                           Delta.gecv = self$Delta.gecv)
    jacknife.gecr.est = function(delta.mat){
      h2.beta.est = sum(delta.mat$Delta.h2.beta);
      h2.gamma.est = sum(delta.mat$Delta.h2.gamma);
      GeCv.est = sum(delta.mat$Delta.gecv);
      return( GeCv.est/sqrt(h2.beta.est *  h2.gamma.est) )
    }
    jacknife.vec = vector()
    for( i in 1:nrow(delta.mat)){
      jacknife.vec = c(jacknife.vec, jacknife.gecr.est(delta.mat[-i,]) )
    }
    
    # ----------------------
    estimator.var = length(Delta)*var(Delta ) + (sum(Delta) -  estimator)^2/sum(!is.na(Delta))
    
    # plugest = self$plug() - mu.y*mu.z
    return(c(estimator, estimator.var))
  }
))
