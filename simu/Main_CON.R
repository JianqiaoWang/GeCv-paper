#system(paste0("plink --bfile ", geno.file ," --make-grm-bin --out ", (geno.file) ))
library(data.table)
library(dplyr)
library(tidyverse)
library(parallel)
library(glmnet)
library(bigstatsr)
library(bigsnpr)

source("GeCr.R")
source("Misc.R")
source("GeCv_bigsnpr.R")

#------------------- method part --------------------------------------
# specify the input parameter
prop.beta <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
prop.gamma <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
GeCovaraince <-as.numeric(commandArgs(trailingOnly = TRUE)[4])
prop.o <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
mixedcoef <- commandArgs(trailingOnly = TRUE)[6]
overlap <- commandArgs(trailingOnly = TRUE)[7]

prop.beta = 0.05 ; prop.gamma = 0.05; GeCovaraince=0.2;prop.o=0.05;mixedcoef ="normal";overlap = "overlap"


i = as.numeric(Sys.getenv("LSB_JOBINDEX"))
# If i = 0 set it to 1, useful for when you're de-bugging code in an interactive session
if(is.na(i)| i ==0){i = 1}

# files directory
#geno.file = "../../Rose_plink_data/rose_all"
if(is.na(mixedcoef)| mixedcoef == "normal" ){

  file.dir = paste0("./Out/L_L/",overlap,"/",
                    prop.beta,"_",prop.gamma,"_",
                    prop.o, "_", GeCovaraince,"/")
}else{

  file.dir = paste0("./Out/L_L/",overlap,"/",
                    prop.beta,"_",prop.gamma,"_",
                    prop.o, "_", GeCovaraince, "_", mixedcoef,"/")
}

print(file.dir)
geno.file = "./plink/rose_all_control"
geno <- bigsnpr::snp_attach(paste0(geno.file, ".rds"))
G <- geno$genotypes

#-------------------------------------
output.dir = paste0(file.dir,i,"/")
pheno.file = paste0(output.dir, "outcome.pheno")
#--------------- function --------------------
SimuBigstats.cv(geno.file = geno.file,
             pheno.file = pheno.file,
             output.dir = output.dir,
             ncores = 1)

#SimuMtg2(geno.file = geno.file,
#             pheno.file = pheno.file,
#             output.dir = output.dir,
#             overlap = overlap)

SimuGCTA(geno.file = geno.file,
         pheno.file = pheno.file,
         output.dir = output.dir,
         overlap = overlap)
