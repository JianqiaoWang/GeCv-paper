library(data.table)
library(dplyr)
library(tidyverse)
library(parallel)
library(glmnet)
source("Misc.R")
source("CoefGen.R")
source("pheno.fun.R")

#--------------Create a vector of seeds for the normal batch job
i = as.numeric(Sys.getenv("LSB_JOBINDEX"))
# If i = 0 set it to 1, useful for when you're de-bugging code in an interactive session
if(is.na(i)| i ==0){i = 1}
#--------------Create a vector of seeds for array job
set.seed(i)
# Parameter input: genotype file, phenotype file and output directory
prop.beta <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
prop.gamma <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
GeCovaraince <-as.numeric(commandArgs(trailingOnly = TRUE)[4])
prop.o <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
mixedcoef <- commandArgs(trailingOnly = TRUE)[6]
overlap <- commandArgs(trailingOnly = TRUE)[7]
prop.beta = 0.02 ; prop.gamma = 0.002; GeCovaraince=0.2;prop.o=0.002;
mixedcoef ="LDAKs3";overlap = "nooverlap"
if(is.na(mixedcoef)| mixedcoef == "normal" ){

  file.dir = paste0("./Out/L_L/",overlap,"/",
                    prop.beta,"_",prop.gamma,"_",
                    prop.o, "_", GeCovaraince,"/")
}else{

  file.dir = paste0("./Out/L_L/", overlap,"/",
                    prop.beta,"_",prop.gamma,"_",
                    prop.o, "_", GeCovaraince, "_", mixedcoef,"/")
}

print(file.dir)
dir.create(file.dir, showWarnings = T, recursive = T)

# parameter settings
h2.Y = 0.5;  h2.Z = 0.5;
Var.Z = 1 ;  Var.Y = 1
sd.ep1 = sqrt(Var.Y * (1 - h2.Y))
sd.ep2 = sqrt(Var.Z * (1 - h2.Z))


#-------------- Function for generating the outcomes -----------------
geno.file = "./plink/rose_all_control"
causal.snp = snpStats::read.plink("./plink/rose_all_control_pruned")
geno = causal.snp$genotypes %>% as("numeric") #%>% dplyr::slice(match())
geno = na.process(geno)
#--------------- Generate the coefficient vectors    -----------------
K = ncol(geno)
d1 = round( (prop.beta - prop.o) * K)
d2 = round( (prop.gamma - prop.o) * K)
d3 = round( prop.o * K)
X <-rep(0, K);
X[sample(1:round(K/3), d1 )] = 1
X[sample(round(2*K/3):K, d2)] = 2
X[sample(which(X ==0), d3 )] = 3
#---------------- Generate the coefficient vectors    -----------------
Coef = Coef.Prop.new(X = X, geno = geno,   h2.Y,
                     sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                     h2.Z, GeCovaraince)

if(mixedcoef == "mixed"){
  Coef = Coef.Prop.new.mixed(X = X, geno = geno,
                             h2.Y,
                             sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                             h2.Z, GeCovaraince)
}

if(mixedcoef == "LDAK"){
  ID = colnames(geno)
  LDSC = data.table::fread("./plink/rose_all_control.l2.ldscore.gz") %>% dplyr::slice(match(ID, SNP))
  Coef = Coef.Prop.new.LDAK(X = X, geno = geno,   h2.Y,
                            sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                            h2.Z, GeCovaraince, LDW = LDSC$L2 )
}

if(mixedcoef == "LDAKs1"){
  ID = colnames(geno)
  LDSC = data.table::fread("./plink/rose_all_control.l2.ldscore.gz") %>% dplyr::slice(match(ID, SNP))
  Coef = Coef.Prop.new.LDAK(X = X, geno = geno, alpha = -0.75,  h2.Y,
                            sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                            h2.Z, GeCovaraince, LDW = LDSC$L2 )
}

if(mixedcoef == "LDAKs2"){
  ID = colnames(geno)
  LDSC = data.table::fread("./plink/rose_all_control.l2.ldscore.gz") %>% dplyr::slice(match(ID, SNP))
  Coef = Coef.Prop.new.LDAK(X = X, geno = geno, alpha = 0.75,  h2.Y,
                            sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                            h2.Z, GeCovaraince, LDW = 1/LDSC$L2 )
}

if(mixedcoef == "LDAKs3"){
  ID = colnames(geno)
  LDSC = data.table::fread("./plink/rose_all_control.l2.ldscore.gz") %>% dplyr::slice(match(ID, SNP))
  Coef = Coef.Prop.new.LDAK(X = X, geno = geno, alpha = -0.75,  h2.Y,
                            sigma3 =  sqrt(prop.beta/prop.o * 0.4),
                            h2.Z, GeCovaraince, LDW = 1/LDSC$L2 )
}

pr_Y = geno %*% Coef$Beta
pr_Z = geno %*% Coef$Gamma1

if(mixedcoef == "DOM" | mixedcoef == "INTER"){
pr_Y = get_pr(geno = geno, coef.vec = Coef$Beta, model = mixedcoef, h2 = h2.Y)
}

if(mixedcoef == "DOM2" ){
  pr_Y = get_pr(geno = geno, coef.vec = Coef$Beta, model = "DOM", h2 = h2.Y)
  pr_Z = get_pr(geno = geno, coef.vec = Coef$Gamma1, model = "DOM", h2 = h2.Z)
}

if(mixedcoef == "INTER2" ){
  pr_Y = get_pr(geno = geno, coef.vec = Coef$Beta, model = "INTER", h2 = h2.Y)
  pr_Z = get_pr(geno = geno, coef.vec = Coef$Gamma1, model = "INTER", h2 = h2.Z)
}

TruePara = cov(pr_Y, pr_Z);
var(pr_Y); var(pr_Z);TruePara;
saveRDS(list(Coef, TruePara), file = paste0(file.dir, "Coef.rds"))
#----------------- Generate the simulated outocmes
simPheno = function(i){
  #----------------------------------------------------------
  # Create the ouput directory
  output.dir = paste0(file.dir,i,"/")
  dir.create(output.dir, showWarnings = T, recursive = T)
  #------------------------Generate phenotype file--------------
  pheno.file = paste0(output.dir, "outcome.pheno")
  mgt.pheno.file = paste0(output.dir, "outcome.dat")


  rho = 0.4
  temp = mvtnorm::rmvnorm(nrow(geno), mean = rep(0,2),
                          sigma = matrix(c(sd.ep1^2, rho*sd.ep1*sd.ep2,
                                           rho*sd.ep1*sd.ep2, sd.ep2^2),2,2))

  Y1 = pr_Y + temp[,1]
  Y2 = pr_Z + temp[,2]
  # ------------- generate the phenotype file ----------------
  IID = rownames(geno)
  I = sample(IID, 8000);
  if(overlap == "overlap") {I.y = I; I.z = I;}
  if(overlap == "nooverlap") {I.y = I[1:4000]; I.z = I[4001:8000];}

  Y1[setdiff( rownames(geno) ,I.y),1] = NA;
  Y2[setdiff(rownames(geno) ,I.z),1] = NA;
  Fam.file = fread(paste0(geno.file, ".fam"));

  split.col = sample(c("val","train", "test"), size = length(Y1),
                     replace = TRUE, prob = c(0.2, 0.7, 0.1) )
  split.col.2 = split.col

  FID = Fam.file$V1 ; IID = Fam.file$V2;

  Outcome.df = data.frame(FID = FID, IID = IID, Y1 = Y1, Y2 = Y2,
                          Pred = rnorm(length(FID)),  split = split.col, split2 = split.col.2)

  Outcome.df %>% dplyr::filter(IID %in% I) %>%
    data.table::fwrite(file = pheno.file, quote = F,
                       sep = "\t", col.names=T, na = "NA")

  Outcome.df %>% dplyr::select(FID, IID, Y1, Y2) %>%
    data.table::fwrite(file = mgt.pheno.file, quote = F,
                       sep = "\t", col.names=F, na = "NA")

}

mclapply(1:200, simPheno, mc.cores = 5)

