get_pr <- function(
  geno,
  coef.vec,
  model,
  h2,
  K = 0.5                                   ## prevalence
){

  if (model == "linear") {
    # only linear
    pr_mean = geno %*% coef.vec

  }
   if (model == "DOM"){

    # dominance
    pr_mean <- (geno != 0) %*% coef.vec * sqrt(2)
    # standardize
    pr_mean <- pr_mean / sd(pr_mean) * sqrt(h2)

  }

  if (model == "INTER"){
    non.zero.ind <- which(coef.vec != 0)
    ind3 = non.zero.ind
    set.seed(234)
    ind3.rand = sample(ind3)
    coef.sign = sign(coef.vec)
    coef.size = abs(coef.vec)
    
    pr1 = sqrt(geno) %*% coef.vec
    pr2 = (geno)^3 %*% coef.vec
    pr3 = (geno[, ind3] * geno[, ind3.rand]  ) %*%  (coef.vec[ind3] * coef.size[ind3.rand])
    pr_mean <- pr1/sd(pr1) + pr2/sd(pr2) + pr3/sd(pr3)
    pr_mean <- pr_mean / sd(pr_mean) * sqrt(h2)
  }

  if (model == "COMP"){
    non.zero.ind <- which(coef.vec != 0)
    sets <- split( non.zero.ind, sample(rep_len(1:3, length(non.zero.ind))))
    # only linear
    ind1 <- sets[[1]]
    pr_mean <- geno[, ind1] %*% coef.vec[ind1]
    # recessive / dominant
    ind2 <- sets[[2]]
    pr_mean <- pr_mean + (geno[, ind2] != 0) %*% coef.vec[ind2]
    # interactions
    ind3 <-sets[[3]]
    pr_mean <- pr_mean + (geno[, ind3[1:(length(ind3) - 1)]] * geno[, ind3[2:(length(ind3))]]) %*%
      coef.vec[ind3[1:(length(ind3) - 1)]]

    # standardize
    pr_mean <- pr_mean / sd(pr_mean) * sqrt(h2)

  }


  if (model == "logit") {
    # only linear
    pr_mean = exp(geno %*% coef.vec)/(1 + exp(geno %*% coef.vec))
  }
 return(pr_mean)
}


get_pr_simple <- function(
  geno,
  coef.vec,
  model,
  h2,
  K = 0.5                                   ## prevalence
){

  if (model == "linear") {
    # only linear
    pr_mean = geno %*% coef.vec

  }
  if (model == "COMP"){
    non.zero.ind <- which(coef.vec != 0)
    sets <- split( non.zero.ind, sample(rep_len(1:2, length(non.zero.ind))))
    # only linear
    ind1 <- sets[[1]]
    pr_mean <- geno[, ind1] %*% coef.vec[ind1]
    # recessive / dominant
    ind2 <- sets[[2]]
    pr_mean <- pr_mean + (geno[, ind2] != 0) %*% coef.vec[ind2] * sqrt(2)
    # interactions
    #ind3 <-sets[[3]]
    #pr_mean <- pr_mean + (geno[, ind3[1:(length(ind3) - 1)]] * geno[, ind3[2:(length(ind3))]]) %*%
    #  coef.vec[ind3[1:(length(ind3) - 1)]]

    # standardize
    pr_mean <- pr_mean / sd(pr_mean) * sqrt(h2)

  }


  if (model == "logit") {
    # only linear
    pr_mean = exp(geno %*% coef.vec)/(1 + exp(geno %*% coef.vec))
  }
 return(pr_mean)
}


expit = function(x){
  exp(x)/(1 + exp(x))
}

linear = function(x){
  return(x)
}

expn = function(x){
  return(x)
}


#
# get_pheno <- function( pr_y, pr_z, outcome.1, outcome.2 ){
#
#
# }
