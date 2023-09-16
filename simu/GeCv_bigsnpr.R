SimuBigstats = function(geno.file,
                        pheno.file,
                        output.dir,
                        ncores = 1){

library(bigstatsr)

Outcome <- read.delim(pheno.file, header=TRUE)

  Y1.train = Outcome$Y1[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train = Outcome$Y2[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  Y1.train.IID = Outcome$IID[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train.IID = Outcome$IID[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  ind.train.Y1 = match(Y1.train.IID, geno$fam$sample.ID)
  ind.train.Y2 = match(Y2.train.IID, geno$fam$sample.ID)
  ind = match(Outcome$IID, geno$fam$sample.ID)
#-----------------------------------------------------------------------
a <- Sys.time()

FUN.1 <- if (all(Y1.train %in% 0:1)) big_spLogReg else big_spLinReg
FUN.2 <- if (all(Y2.train %in% 0:1)) big_spLogReg else big_spLinReg

fit1 = FUN.1(
  X = G,
  Y1.train,
  ind.train = ind.train.Y1,
  alphas = 1,
  power_scale = 1,
  power_adaptive = 0,
  K = 10,
  ncores = ncores)

Y.hat = predict(fit1, X = G, ind.row = ind)

fit2 = FUN.2(
  X = G,
  Y2.train,
  ind.train = ind.train.Y2,
  alphas = 1,
  power_scale = 1,
  power_adaptive = 0,
  K = 10,
  ncores = ncores)
b <- Sys.time()
msg1 = (paste0("training time is", b-a))
print(msg1)
#---------------------------------------------------------------------------
Z.hat = predict(fit2, X = G, ind.row = ind)
Dataset = data.frame(Y = Outcome$Y1,
                     Z = Outcome$Y2,
                     Y.hat = Y.hat,
                     Z.hat = Z.hat )
Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
Dataset$eps.hat[is.na(Dataset$Y)] = 0; Dataset$v.hat[is.na(Dataset$Z)] = 0;
n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
attach(Dataset)
GePa = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat, v.hat =  Dataset$v.hat)
GePa.full = GePa
result.full = c(GePa$h2.beta.est/var(Dataset$Y, na.rm = T),
           GePa$h2.beta.est.se/var(Dataset$Y, na.rm = T),
           GePa$h2.gamma.est/var(Dataset$Z, na.rm = T),
           GePa$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
           GePa$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
           GePa$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
detach(Dataset)


ind.split = Outcome$split == "test"|Outcome$split2 == "test"
Dataset = Dataset[ind.split, ]
n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
attach(Dataset)
GePa.split = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat = Dataset$Y.hat,
                Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat, v.hat =  Dataset$v.hat)
result.split = c(GePa.split$h2.beta.est/var(Dataset$Y, na.rm = T),
                 GePa.split$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                 GePa.split$h2.gamma.est/var(Dataset$Z, na.rm = T),
                 GePa.split$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                 GePa.split$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                 GePa.split$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
detach(Dataset)
saveRDS(list(GePa.full, result.full, GePa.split, result.split, Dataset),
        file = paste0(output.dir, "bigstat.GeCv.rds" ))
return(msg1)
}


SimuBigstats.cv = function(geno.file,
                        pheno.file,
                        output.dir,
                        ncores = 1,
                        scale_var = T){

Outcome <- read.delim(pheno.file, header=TRUE)

  Y1.train = Outcome$Y1[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train = Outcome$Y2[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  Y1.train.IID = Outcome$IID[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train.IID = Outcome$IID[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  ind.train.Y1 = match(Y1.train.IID, geno$fam$sample.ID)
  ind.train.Y2 = match(Y2.train.IID, geno$fam$sample.ID)
  ind = match(Outcome$IID, geno$fam$sample.ID)
#-----------------------------------------------------------------------
a <- Sys.time()

FUN.1 <- if (all(Y1.train %in% 0:1)) big_spLogReg.cv else big_spLinReg.cv
FUN.2 <- if (all(Y2.train %in% 0:1)) big_spLogReg.cv else big_spLinReg.cv

fit1 = FUN.1(
  X = G,
  Y1.train,
  ind.train = ind.train.Y1,
  alphas = 1,
  power_scale = 1,
  power_adaptive = 0,
  K = 10,
  warn = FALSE,
  ncores = ncores)

Y.hat = predict.cv(fit1, X = G, ind.row = ind)

fit2 = FUN.2(
  X = G,
  Y2.train,
  ind.train = ind.train.Y2,
  alphas = 1,
  power_scale = 1,
  power_adaptive = 0,
   warn = FALSE,
  K = 10,
  ncores = ncores)

b <- Sys.time()
msg1 = (paste0("training time is", b-a))
print(msg1)
#---------------------------------------------------------------------------
Z.hat = predict.cv(fit2, X = G, ind.row = ind)
Dataset = data.frame(Y = Outcome$Y1,
                     Z = Outcome$Y2,
                     Y.hat = Y.hat,
                     Z.hat = Z.hat )
Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
Dataset$eps.hat[is.na(Dataset$Y)] = 0; Dataset$v.hat[is.na(Dataset$Z)] = 0;
n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
attach(Dataset)
GePa = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat, v.hat =  Dataset$v.hat)
GePa.full = GePa
if(scale_var){
result.full = c(GePa$h2.beta.est/var(Dataset$Y, na.rm = T),
           GePa$h2.beta.est.se/var(Dataset$Y, na.rm = T),
           GePa$h2.gamma.est/var(Dataset$Z, na.rm = T),
           GePa$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
           GePa$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
           GePa$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
}else{
  
  result.full = c(GePa$h2.beta.est,
                  GePa$h2.beta.est.se,
                  GePa$h2.gamma.est,
                  GePa$h2.gamma.est.se,
                  GePa$GeCv.est,
                  GePa$GeCv.est.se)
} 
detach(Dataset)


ind.split = Outcome$split == "test"|Outcome$split2 == "test"
Dataset = Dataset[ind.split, ]
n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
attach(Dataset)
GePa.split = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat = Dataset$Y.hat,
                Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat, v.hat =  Dataset$v.hat)
if(scale_var){
result.split = c(GePa.split$h2.beta.est/var(Dataset$Y, na.rm = T),
                 GePa.split$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                 GePa.split$h2.gamma.est/var(Dataset$Z, na.rm = T),
                 GePa.split$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                 GePa.split$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                 GePa.split$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
}else{
  result.split = c(GePa.split$h2.beta.est,
                   GePa.split$h2.beta.est.se,
                   GePa.split$h2.gamma.est,
                   GePa.split$h2.gamma.est.se,
                   GePa.split$GeCv.est,
                   GePa.split$GeCv.est.se)
  
}

attr(fit1, "pf") <- NULL
attr(fit1, "ind.col") <- NULL
attr(fit1, "ind.sets") <- NULL
attr(fit1, "base") <- NULL

detach(Dataset)
saveRDS(list(GePa.full, result.full, GePa.split, result.split, Dataset, fit1),
        file = paste0(output.dir, "bigstat.cv.GeCv.rds" ))
return(msg1)
}


CCSimuBigstats = function(geno.file,
                        pheno.file,
                        output.dir,
                        prev.Y = NULL,
                        prev.Z = NULL,
                        ncores = 1){

  Outcome <- read.delim(pheno.file, header=TRUE)

  Y1.train = Outcome$Y1[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train = Outcome$Y2[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  Y1.train.IID = Outcome$IID[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train.IID = Outcome$IID[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  ind.train.Y1 = match(Y1.train.IID, geno$fam$sample.ID)
  ind.train.Y2 = match(Y2.train.IID, geno$fam$sample.ID)
  ind = match(Outcome$IID, geno$fam$sample.ID)
  #-----------------------------------------------------------------------
  a <- Sys.time()

  FUN.1 <- if (all(Y1.train %in% 0:1)) big_spLogReg else big_spLinReg
  FUN.2 <- if (all(Y2.train %in% 0:1)) big_spLogReg else big_spLinReg

  fit1 = FUN.1(
    X = G,
    Y1.train,
    ind.train = ind.train.Y1,
    alphas = 1,
    power_scale = 1,
    power_adaptive = 0,
    K = 10,
    ncores = ncores)

  fit.summary = summary(fit1)
  intercept = fit.summary$intercept
  if(!is.null(prev.Y)){
    p1 = prev.Y; p0 = 1 - prev.Y;
    pi1 = sum(Y1.train == 1, na.rm = T)/ sum(!is.na(Y1.train));
    pi0 = sum(Y1.train == 0, na.rm = T)/ sum(!is.na(Y1.train ));
    w0.Y = p0/pi0
    w1.Y = p1/pi1
    intercept = intercept + log(w1.Y/w0.Y)
  }
  beta = unlist(fit.summary$beta)
  ind.nozero <- which(beta != 0); ind.col = attr(fit1, "ind.col")
  scores <- big_prodVec(X = G, beta[ind.nozero], ind.row = ind,
                        ind.col = ind.col[ind.nozero],
                        ncores = ncores)
  Y.hat = expit( intercept + scores )


  #Y.hat = predict(fit1, X = G, ind.row = ind)


  fit2 = FUN.2(
    X = G,
    Y2.train,
    ind.train = ind.train.Y2,
    alphas = 1,
    power_scale = 1,
    power_adaptive = 0,
    K = 10,
    ncores = ncores)
  b <- Sys.time()
  msg1 = (paste0("training time is", b-a))
  print(msg1)

  fit.summary = summary(fit2)
  intercept = fit.summary$intercept
  if(!is.null(prev.Z)){
    p1 = prev.Z; p0 = 1 - prev.Z;
    pi1 = sum(Y2.train == 1, na.rm = T)/ sum(!is.na(Y2.train));
    pi0 = sum(Y2.train == 0, na.rm = T)/ sum(!is.na(Y2.train ));
    w0.Z = p0/pi0
    w1.Z = p1/pi1
    intercept = intercept + log(w1.Z/w0.Z)
  }
  beta = unlist(fit.summary$beta)
  ind.nozero <- which(beta != 0) ; ind.col = attr(fit2, "ind.col")
  scores <- big_prodVec(X = G, beta[ind.nozero], ind.row = ind,
                        ind.col = ind.col[ind.nozero], ncores = ncores)
  Z.hat = expit( intercept + scores )
  #---------------------------------------------------------------------------
  #Z.hat = predict(fit2, X = G,
  #                ind.row = ind)


  Dataset = data.frame(Y = Outcome$Y1,
                       Z = Outcome$Y2,
                       Y.hat = Y.hat,
                       Z.hat = Z.hat )
  Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
  Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
  Dataset$eps.hat[is.na(Dataset$Y)] = 0;
  Dataset$v.hat[is.na(Dataset$Z)] = 0;
  Dataset$weight = 1

  Dataset$weight[which(Dataset$Y == 1)] = w1.Y
  Dataset$weight[which(Dataset$Y == 0)] = w0.Y
  Dataset$weight[which(Dataset$Z == 1)] = w1.Z
  Dataset$weight[which(Dataset$Z == 0)] = w0.Z

  n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
  attach(Dataset)
  GePa = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                  Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat,
                  v.hat =  Dataset$v.hat, weight = Dataset$weight )
  GePa.full = GePa
  result.full = c(GePa$h2.beta.est/var(Dataset$Y, na.rm = T),
                  GePa$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                  GePa$h2.gamma.est/var(Dataset$Z, na.rm = T),
                  GePa$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                  GePa$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                  GePa$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
  detach(Dataset)


  ind.split = Outcome$split == "test"|Outcome$split2 == "test"
  Dataset = Dataset[ind.split, ]
  n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
  attach(Dataset)
  GePa.split = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                        Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat,
                        v.hat =  Dataset$v.hat,
                        weight = Dataset$weight )
  result.split = c(GePa.split$h2.beta.est/var(Dataset$Y, na.rm = T),
                   GePa.split$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                   GePa.split$h2.gamma.est/var(Dataset$Z, na.rm = T),
                   GePa.split$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                   GePa.split$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                   GePa.split$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
  detach(Dataset)
  saveRDS(list(GePa.full, result.full, GePa.split, result.split, Dataset),
          file = paste0(output.dir, "bigstat.GeCv.rds" ))
  return(msg1)
}




CCSimuBigstats.True = function(geno.file,
                        pheno.file,
                        output.dir,
                        prev.Y = NULL,
                        prev.Z = NULL,
                        Coef = NULL,
                        ncores = 1){

  Outcome <- read.delim(pheno.file, header=TRUE)

  Y1.train = Outcome$Y1[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train = Outcome$Y2[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  Y1.train.IID = Outcome$IID[ (!is.na(Outcome$Y1)) & Outcome$split %in%c("train", "val") ]
  Y2.train.IID = Outcome$IID[ (!is.na(Outcome$Y2)) & Outcome$split2 %in%c("train", "val") ]
  ind.train.Y1 = match(Y1.train.IID, geno$fam$sample.ID)
  ind.train.Y2 = match(Y2.train.IID, geno$fam$sample.ID)
  ind = match(Outcome$IID, geno$fam$sample.ID)
  #-----------------------------------------------------------------------
  a <- Sys.time()

  FUN.1 <- if (all(Y1.train %in% 0:1)) big_spLogReg else big_spLinReg
  FUN.2 <- if (all(Y2.train %in% 0:1)) big_spLogReg else big_spLinReg

  fit1 = FUN.1(
    X = G,
    Y1.train,
    ind.train = ind.train.Y1,
    alphas = 1,
    power_scale = 1,
    power_adaptive = 0,
    K = 10,
    ncores = ncores)

  fit.summary = summary(fit1)
  intercept = fit.summary$intercept
  if(!is.null(prev.Y)){
    p1 = prev.Y; p0 = 1 - prev.Y;
    pi1 = sum(Y1.train == 1, na.rm = T)/ sum(!is.na(Y1.train));
    pi0 = sum(Y1.train == 0, na.rm = T)/ sum(!is.na(Y1.train ));
    w0.Y = p0/pi0
    w1.Y = p1/pi1
    intercept = intercept + log(w1.Y/w0.Y)
  }
  beta = unlist(fit.summary$beta)
  beta = Coef$Beta
  ind.nozero <- which(beta != 0); ind.col = attr(fit1, "ind.col")
  scores <- big_prodVec(X = G, beta[ind.nozero], ind.row = ind,
                        ind.col = ind.col[ind.nozero],
                        ncores = ncores)
  Y.hat = expit( intercept + scores )


  #Y.hat = predict(fit1, X = G, ind.row = ind)


  fit2 = FUN.2(
    X = G,
    Y2.train,
    ind.train = ind.train.Y2,
    alphas = 1,
    power_scale = 1,
    power_adaptive = 0,
    K = 10,
    ncores = ncores)
  b <- Sys.time()
  msg1 = (paste0("training time is", b-a))
  print(msg1)

  fit.summary = summary(fit2)
  intercept = fit.summary$intercept
  if(!is.null(prev.Z)){
    p1 = prev.Z; p0 = 1 - prev.Z;
    pi1 = sum(Y2.train == 1, na.rm = T)/ sum(!is.na(Y2.train));
    pi0 = sum(Y2.train == 0, na.rm = T)/ sum(!is.na(Y2.train ));
    w0.Z = p0/pi0
    w1.Z = p1/pi1
    intercept = intercept + log(w1.Z/w0.Z)
  }
  beta = unlist(fit.summary$beta)
  beta = Coef$Gamma1
  ind.nozero <- which(beta != 0) ; ind.col = attr(fit2, "ind.col")
  scores <- big_prodVec(X = G, beta[ind.nozero], ind.row = ind,
                        ind.col = ind.col[ind.nozero], ncores = ncores)
  Z.hat = expit( intercept + scores )
  #---------------------------------------------------------------------------
  #Z.hat = predict(fit2, X = G,
  #                ind.row = ind)


  Dataset = data.frame(Y = Outcome$Y1,
                       Z = Outcome$Y2,
                       Y.hat = Y.hat,
                       Z.hat = Z.hat )
  Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
  Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
  Dataset$eps.hat[is.na(Dataset$Y)] = 0;
  Dataset$v.hat[is.na(Dataset$Z)] = 0;
  Dataset$weight = 1

  Dataset$weight[which(Dataset$Y == 1)] = w1.Y
  Dataset$weight[which(Dataset$Y == 0)] = w0.Y
  Dataset$weight[which(Dataset$Z == 1)] = w1.Z
  Dataset$weight[which(Dataset$Z == 0)] = w0.Z

  n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
  attach(Dataset)
  GePa = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                  Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat,
                  v.hat =  Dataset$v.hat, weight = Dataset$weight )
  GePa.full = GePa
  result.full = c(GePa$h2.beta.est/var(Dataset$Y, na.rm = T),
                  GePa$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                  GePa$h2.gamma.est/var(Dataset$Z, na.rm = T),
                  GePa$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                  GePa$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                  GePa$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
  detach(Dataset)


  ind.split = Outcome$split == "test"|Outcome$split2 == "test"
  Dataset = Dataset[ind.split, ]
  n.z = sum(!is.na(Dataset$Z)); n.y = sum(!is.na(Dataset$Y));
  attach(Dataset)
  GePa.split = GeCr$new(N.Y = n.y, N.Z = n.z, Y.hat =  Dataset$Y.hat,
                        Z.hat =  Dataset$Z.hat, eps.hat =  Dataset$eps.hat,
                        v.hat =  Dataset$v.hat,
                        weight = Dataset$weight )
  result.split = c(GePa.split$h2.beta.est/var(Dataset$Y, na.rm = T),
                   GePa.split$h2.beta.est.se/var(Dataset$Y, na.rm = T),
                   GePa.split$h2.gamma.est/var(Dataset$Z, na.rm = T),
                   GePa.split$h2.gamma.est.se/var(Dataset$Z, na.rm = T),
                   GePa.split$GeCv.est/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)),
                   GePa.split$GeCv.est.se/sqrt(var(Dataset$Z, na.rm = T)*var(Dataset$Y, na.rm = T)))
  detach(Dataset)
  saveRDS(list(GePa.full, result.full, GePa.split, result.split, Dataset),
          file = paste0(output.dir, "bigstat.GeCv.rds" ))
  return(msg1)
}



SimubigLasso = function(geno.file,
                        pheno.file,
                        output.dir,
                        ncores = 1){

  #------------------------------------------------------------------------------
  Outcome <- read.delim(pheno.file, header=TRUE)
  Y1 = Outcome$Y1; Y2 = Outcome$Y2
  ind.train.Y1 = which(!is.na(Y1) & Outcome$split %in%c("train", "val")  )
  ind.train.Y2 = which(!is.na(Y2) & Outcome$split %in%c("train", "val"))
  bim = fread(paste0(geno.file, ".bim"))
  #------------------------------------------------------------------------------
  a <- Sys.time()
  G <- bigmemory::attach.big.matrix(paste0(geno.file, ".desc"))
  family.1 <- if (all(Y1[ind.train.Y1] %in% 0:1)) "binomial" else "gaussian"
  family.2 <- if (all(Y2[ind.train.Y2] %in% 0:1)) "binomial" else "gaussian"

  fit1 = cv.biglasso2(X = G, Y1, row.idx  = ind.train.Y1, ncores = ncores,  nfolds = 10,
                     family = family.1, seed = 123)


  Y.hat = predict(fit1, X = G,type = "response")

  fit2 = cv.biglasso2(X = G, Y2, row.idx  = ind.train.Y2, ncores = ncores, family = family.2,
                       nfolds = 10, seed = 123 )
  #fit2 = cv.biglasso(X = G, Y2,  ncores = ncores, family = family.2,
  #                   nfolds = 10, seed = 123 )
  Z.hat = predict(fit2, X = G, type = "response")
  #-----------------------------------------------------------------------------
  b <- Sys.time();  print(paste0("training time is", b-a))
  #-----------------------------------------------------------------------------
  Dataset = data.frame(Y = Outcome$Y1,  Z = Outcome$Y2,   Y.hat = as.matrix(Y.hat),  Z.hat = as.matrix(Z.hat) )
  Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
  Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
  Dataset$eps.hat[is.na(Dataset$Y)] = 0; Dataset$v.hat[is.na(Dataset$Z)] = 0;
  ind.split = Outcome$split == "test"|Outcome$split2 == "test"
  res.full = res.fun(Dataset[-ind.split, ])
  res.split = res.fun(Dataset[ind.split, ])
  #--------------------- POST OLS -----------------------------------------------
  beta.hat = coef(fit1, s= "lambda.min")
  gamma.hat = coef(fit2, s= "lambda.min")
  if(any(beta.hat@i != 0) ){
    G.sub.Y =   G[, beta.hat@i] %>% as.data.frame
  }else{
    G.sub.Y = as.data.frame(rep(1, length(Y1)))
  }
  lm1 = lm(Y1 ~ as.matrix(G.sub.Y), subset = ind.train.Y1)

  if(any(gamma.hat@i != 0)){
    G.sub.Z =   G[, gamma.hat@i] %>% as.data.frame
  }else{
    G.sub.Z =   as.data.frame(rep(1, length(Y2)))
  }
  lm2 = lm(Y2 ~ as.matrix(G.sub.Z), subset = ind.train.Y2 )

  Y.hat.OLS = predict(lm1, newdata = G.sub.Y,type = "response")
  Z.hat.OLS = predict(lm2, newdata = G.sub.Z, type = "response")
  Dataset = data.frame(Y = Outcome$Y1,  Z = Outcome$Y2,
                       Y.hat = Y.hat.OLS,  Z.hat = Z.hat.OLS)
  Dataset$eps.hat <- Dataset$Y - Dataset$Y.hat;
  Dataset$v.hat <- Dataset$Z - Dataset$Z.hat;
  Dataset$eps.hat[is.na(Dataset$Y)] = 0;
  Dataset$v.hat[is.na(Dataset$Z)] = 0;
  res.full.postOLS = res.fun(Dataset[-ind.split, ])
  res.split.postOLS = res.fun(Dataset[ind.split, ])
  #-----------------

  saveRDS(list(res.full, res.split, res.full.postOLS, res.split.postOLS),
          file = paste0(output.dir, "biglasso.split.GeCv.rds" ))
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

predict.modify = function(fit, X , ind.row){
  fit.summary = summary(fit)
  intercept = fit.summary$intercept
  beta = unlist(fit.summary$beta)
  ind.nozero <- which(beta != 0) ; ind.col = attr(fit, "ind.col")
  scores <- big_prodVec(X = G, beta[ind.nozero], ind.row = ind,
                        ind.col = ind.col[ind.nozero], ncores = ncores)

}

