#' @title Generate coefficients
#' @description Generate coefficients



Coef.Prop.new = function(X, geno, h2.Y, h2.Z,
                       GeCovaraince,
                       mu3 = 0, sigma3 = 1, alpha = 0, Var.Y = 1, Var.Z = 1){

  p = ncol(geno)

  beta.coef = rep(0, p);
  beta.coef[X==1] = rnorm( sum(X==1) , 0, 1)
  beta.coef[X==3] = rnorm( sum(X==3) , mu3, sigma3)

  #----------- Generate Beta coef --------------
  beta.coef <- beta.coef*sqrt(h2.Y*Var.Y)/sqrt( var(geno %*% beta.coef ) ) %>% as.vector()
  #----------- Generate eta and gamma coef --------------
  # construct the orthorgonal vectors to the beta
  x.1 <- rep(0, p)
  x.2 <- rep(0, p)
  x.1[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3))
  x.2[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3))
  x.1[X == 3] <- rnorm(sum(X == 3), mu3, sigma3)
  x.2[X == 3] <- rnorm(sum(X == 3), mu3, sigma3)
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.1)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.2)
  ratio1 = -vpd2/vpd1
  beta.ortho.1 <- ( x.1 * c(ratio1) + x.2)/c(1 + ratio1^2)
  #beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1


  x.3 <- rep(0, p)
  x.4 <- rep(0, p)
  x.3[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3))
  x.4[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3))
  x.3[X == 3] <- rnorm(sum(X == 3), mu3, sigma3)
  x.4[X == 3] <- rnorm(sum(X == 3), mu3, sigma3)
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.3)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.4)
  ratio1 = -vpd2/vpd1
  beta.ortho.2 <- ( x.3 * c(ratio1) + x.4)/c(1 + ratio1^2)
  #beta.ortho.2 = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1

  beta.2 = rep(0, p)
  beta.2[X == 2|X == 3] = beta.coef[X == 2|X == 3] # normalize V1
  pd1 =  cov(geno %*% beta.2, geno %*% beta.ortho.1)
  pd2 =  cov(geno %*% beta.2, geno %*% beta.ortho.2)
  a1 = pd2
  a2 = -pd1
  beta.ortho =  beta.ortho.1 * c(a1)  + beta.ortho.2 * c(a2)


  if(GeCovaraince == 0){
    Gamma1 = beta.ortho.1
    Gamma1 <- Gamma1*sqrt(h2.Z*Var.Z)/sqrt( var(geno %*% Gamma1 ) ) %>% as.vector()
  }else{
  #---------------------------------------------------
  beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1
  h2.Y.part = var(geno %*% beta.2)
  a1 <- (sqrt(Var.Y*Var.Z)*GeCovaraince/ (h2.Y.part*Var.Y) )
  a2 <- sqrt(1 - (Var.Y*h2.Y.part)/(Var.Z*h2.Z) * a1^2  )
  Gamma1 <- as.vector(  beta.2 * c(a1) + beta.ortho * c(a2) )
  }

  return(data.frame(Beta = beta.coef, Gamma1 = Gamma1))
}



Coef.Prop.new.mixed = function(X, geno, h2.Y, h2.Z,
                         GeCovaraince,
                         mu3 = 0, sigma3 = 1, alpha = 0, Var.Y = 1, Var.Z = 1){

  p = ncol(geno)

  rmixed = function(n, mu = 0, sigma = 1){
    temp = cbind(rnorm(n), rnorm(n,1,0.5), rnorm(n,-1,1), rbeta(n,1,1))
    id = sample(1:4,n,rep = T,prob = c(.25,.25,.25, 0.25))
    id = cbind(1:n,id)
    return(temp[id])
  }

  beta.coef = rmixed(p);
  beta.coef[X==2] = 0
  beta.coef[X==3] = rmixed(sum(X==3) , mu3, 1)

  #----------- Generate Beta coef --------------
  beta.coef <- beta.coef*sqrt(h2.Y*Var.Y)/sqrt( var(geno %*% beta.coef ) ) %>% as.vector()
  #gamma.coef <- gamma.coef*sqrt(h2.Z*Var.Z)/sqrt( var(geno %*% gamma.coef ) ) %>% as.vector()

  #----------- Generate eta and gamma coef --------------
  # construct the orthorgonal vectors to the beta
  x.1 <- rep(0, p)
  x.2 <- rep(0, p)
  x.1[X == 2|X == 3] <- rmixed(sum(X == 2|X == 3))
  x.2[X == 2|X == 3] <- rmixed(sum(X == 2|X == 3))
  x.1[X == 3] <- rmixed(sum(X == 3), mu3, 1)
  x.2[X == 3] <- rmixed(sum(X == 3), mu3, 1)
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.1)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.2)
  ratio1 = -vpd2/vpd1
  beta.ortho.1 <- ( x.1 * c(ratio1) + x.2)/c(1 + ratio1^2)
  #beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1


  x.3 <- rep(0, p)
  x.4 <- rep(0, p)
  x.3[X == 2|X == 3] <- rmixed(sum(X == 2|X == 3))
  x.4[X == 2|X == 3] <- rmixed(sum(X == 2|X == 3))
  x.3[X == 3] <- rmixed(sum(X == 3), mu3, 1)
  x.4[X == 3] <- rmixed(sum(X == 3), mu3, 1)
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.3)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.4)
  ratio1 = -vpd2/vpd1
  beta.ortho.2 <- ( x.3 * c(ratio1) + x.4)/c(1 + ratio1^2)
  #beta.ortho.2 = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1

  beta.2 = rep(0, p)
  beta.2[X == 2|X == 3] = beta.coef[X == 2|X == 3] # normalize V1
  pd1 =  cov(geno %*% beta.2, geno %*% beta.ortho.1)
  pd2 =  cov(geno %*% beta.2, geno %*% beta.ortho.2)
  a1 = pd2
  a2 = -pd1
  beta.ortho =  beta.ortho.1 * c(a1)  + beta.ortho.2 * c(a2)


  if(GeCovaraince == 0){
    Gamma1 = beta.ortho.1
    Gamma1 <- Gamma1*sqrt(h2.Z*Var.Z)/sqrt( var(geno %*% Gamma1 ) ) %>% as.vector()
  }else{
    #---------------------------------------------------
    beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1
    h2.Y.part = var(geno %*% beta.2)
    a1 <- (sqrt(Var.Y*Var.Z)*GeCovaraince/ (h2.Y.part*Var.Y) )
    a2 <- sqrt(1 - (Var.Y*h2.Y.part)/(Var.Z*h2.Z) * a1^2  )
    Gamma1 <- as.vector(  beta.2 * c(a1) + beta.ortho * c(a2) )
  }
  return(data.frame(Beta = beta.coef, Gamma1 = Gamma1))
}


Coef.Prop.new.LDAK = function(X, geno, h2.Y, h2.Z, GeCovaraince,
                               mu3 = 0, sigma3 = 1, alpha = 0.75, Var.Y = 1, Var.Z = 1,LDW){

  p = ncol(geno)
  MAF = apply(geno, 2, function(x){mean(x, na.rm = T)/2} )
  #-------------------------------

  weight.sd = (MAF *(1 - MAF))^(-alpha) * 1/LDW
  beta.coef = rep(0, p);
  beta.coef[X==1] = rnorm( sum(X==1) , 0, 1) * weight.sd[X==1]
  beta.coef[X==3] = rnorm( sum(X==3) , mu3, sigma3) * weight.sd[X==3]

  #----------- Generate Beta coef --------------
  beta.coef <- beta.coef*sqrt(h2.Y*Var.Y)/sqrt( var(geno %*% beta.coef ) ) %>% as.vector()
  #----------- Generate eta and gamma coef --------------
  # construct the orthorgonal vectors to the beta
  x.1 <- rep(0, p)
  x.2 <- rep(0, p)
  x.1[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3)) * weight.sd[X == 2|X == 3]
  x.2[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3))* weight.sd[X == 2|X == 3]
  x.1[X == 3] <- rnorm(sum(X == 3), mu3, sigma3) * weight.sd[X == 3]
  x.2[X == 3] <- rnorm(sum(X == 3), mu3, sigma3) * weight.sd[X == 3]
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.1)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.2)
  ratio1 = -vpd2/vpd1
  beta.ortho.1 <- ( x.1 * c(ratio1) + x.2)/c(1 + ratio1^2)
  #beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1


  x.3 <- rep(0, p)
  x.4 <- rep(0, p)
  x.3[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3)) * weight.sd[X == 2|X == 3]
  x.4[X == 2|X == 3] <- rnorm(sum(X == 2|X == 3)) * weight.sd[X == 2|X == 3]
  x.3[X == 3] <- rnorm(sum(X == 3), mu3, sigma3) * weight.sd[X == 3]
  x.4[X == 3] <- rnorm(sum(X == 3), mu3, sigma3) * weight.sd[X == 3]
  vpd1 <- cov(geno %*% beta.coef, geno %*% x.3)
  vpd2 <- cov(geno %*% beta.coef, geno %*% x.4)
  ratio1 = -vpd2/vpd1
  beta.ortho.2 <- ( x.3 * c(ratio1) + x.4)/c(1 + ratio1^2)
  #beta.ortho.2 = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1

  beta.2 = rep(0, p)
  beta.2[X == 2|X == 3] = beta.coef[X == 2|X == 3] # normalize V1
  pd1 =  cov(geno %*% beta.2, geno %*% beta.ortho.1)
  pd2 =  cov(geno %*% beta.2, geno %*% beta.ortho.2)
  a1 = pd2
  a2 = -pd1
  beta.ortho =  beta.ortho.1 * c(a1)  + beta.ortho.2 * c(a2)


  if(GeCovaraince == 0){
    Gamma1 = beta.ortho.1
    Gamma1 <- Gamma1*sqrt(h2.Z*Var.Z)/sqrt( var(geno %*% Gamma1 ) ) %>% as.vector()
  }else{
    #---------------------------------------------------
    beta.ortho = beta.ortho*sqrt(Var.Z*h2.Z)/as.numeric(sqrt(var(geno %*% beta.ortho))) # normalize V1
    h2.Y.part = var(geno %*% beta.2)
    a1 <- (sqrt(Var.Y*Var.Z)*GeCovaraince/ (h2.Y.part*Var.Y) )
    if((Var.Y*h2.Y.part)/(Var.Z*h2.Z) * a1^2 < 1){
    a2 <- sqrt(1 - (Var.Y*h2.Y.part)/(Var.Z*h2.Z) * a1^2  )
    }else{
      a2 = 0
    }
    Gamma1 <- as.vector(  beta.2 * c(a1) + beta.ortho * c(a2) )
    Gamma1 <- Gamma1*sqrt(h2.Z*Var.Z)/sqrt( var(geno %*% Gamma1 ) ) %>% as.vector()

  }

  return(data.frame(Beta = beta.coef, Gamma1 = Gamma1))
}

na.process <- function(X){
  apply(X,2,function(x)
  {
    x[is.na(x)] <- mean(x,na.rm=TRUE);
    return(x);
  });
} # Imputation of NA values


bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

