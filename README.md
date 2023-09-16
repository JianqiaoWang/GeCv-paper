Reproducibility Materials
================

This GitHub repository includes materials to reproduce analyses, visualizations, and
tables for the paper "A Regression-based Approach to Robust Estimation and
Inference for Genetic Covariance".

demo.Rmd (important)
  - This document gives an example of how we conduct the simulation based on the proposed method and how to adjust the confidence interval estimations

Simu/
  - contains code that I used for conduting GWAS simulation
  
  Simu/CoefGen.R
  - a function to generate genetic effects
  - prop parameter control for non-zero numbers: sparse or polygeneic
  - could choose normal or LDAK distriubtion model
  - specifcy the distribution of overlapping effects   

  Simu/GeCr.R
  - define the GeCr class
  - Calculate genetic covariance based on the estimated Y.hat, Z.hat, eps.hat, v.hat 
  
  Simu/GeCv_bigsnpr.R + Misc.R
  - function I used to implement penalized regression with bigstatsr and call GCTA

  Simu/Pheno_CON.R
  -  Mian phenotype generation codes based on the pre-specified parameters. Generate phenotype file. 

  Simu/pheno.fun.R
  -  linear component (or "liability") generation functions: "Linear", "DOM", "INTER"
  
  Simu/Main_CON.R
  - for the give genotype file and phenotype file, main codes for simulation.

