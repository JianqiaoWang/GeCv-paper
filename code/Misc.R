
WritePheno = function( Y1, Y2, FID, IID, output.file  ){

  pheno = data.frame(FID = FID, IID = IID,
                     Y1 = Y1, Y2 = Y2)
  fwrite(pheno, file = output.file,
         quote = F, sep = "\t", col.names=T, na = "NA")

}

SimuGCTA = function(geno.file, pheno.file, output.dir){

  output.file = paste0(output.dir, "GCTA.result")

  system(paste0("./gcta64 --reml-bivar --grm ",
                geno.file," --pheno ", pheno.file,
                " --thread-num 20 --out ", output.file ))
  # if(file.exists(paste0(output.file,".hsq"))){
  #   hsq = fread(paste0(output.file,".hsq"), fill = T)
  #   pval = as.numeric(sapply(strsplit(as.character(hsq[hsq$Source == "Pval",2]), "\\("), "[[", 1))
  # }else{ pval = NA }
  # return(pval)
}

# for case-control data

SimuPCGC = function(geno.file, pheno.file, output.file){
  #system(paste0("plink --bfile ", geno.file ," --make-grm-bin --out ", paste0(i,"simu")))

#outputfile="/home/wangjq/Linear-based/GWAS_Simu/Out/B_B/logistic/0.1_0.1_0.1_0.2/2/s1.pheno"

#python S-PCGC/pcgc_sumstats_creator.py \
#--bfile ${genofile} \
#--pheno ${outputfile} \
#--frqfile ${genofile}. \
#--annot ${genofile}. \
#--sync ${genofile}. \
#--prev 0.1 \
#--out /home/wangjq/Linear-based/GWAS_Simu/Out/B_B/logistic/0.1_0.1_0.1_0.2/2/s1

#outputfile="/home/wangjq/Linear-based/GWAS_Simu/Out/B_B/logistic/0.1_0.1_0.1_0.2/2/s2.pheno"

#python S-PCGC/pcgc_sumstats_creator.py \
#--bfile ${genofile} \
#--pheno ${outputfile} \
#--frqfile ${genofile}. \
#--annot ${genofile}. \
#--sync ${genofile}. \
#--prev 0.1 \
#--out /home/wangjq/Linear-based/GWAS_Simu/Out/B_B/logistic/0.1_0.1_0.1_0.2/2/s2

#filedir='/home/wangjq/Linear-based/GWAS_Simu/Out/B_B/logistic/0.1_0.1_0.1_0.2/2/'
#genofile="./Case_Control_plink/geno"

#python S-PCGC/pcgc_main.py \
# --M 13649 \
# --prodr2 ${genofile}. \
# --sumstats ${filedir}s1.\
# --out ${filedir}results

  system(paste0("gcta64 --reml-bivar --reml-bivar-nocove --reml-maxit 100 --grm ",
                geno.file," --pheno ", pheno.file,
                " --reml-bivar-lrt-rg 0 --thread-num 20 --out ", output.file ))
  if(file.exists(paste0(output.file,".hsq"))){
    hsq = fread(paste0(output.file,".hsq"), fill = T)
    pval = as.numeric(sapply(strsplit(as.character(hsq[hsq$Source == "Pval",2]), "\\("), "[[", 1))
  }else{ pval = NA }
  #./mtg2 -p simu.fam -pheno simu.phen --grm simu -out simu  -mod 2
  #test.out -mod 5
  return(pval)
}

SimuMtg2 = function(geno.file, pheno.file, output.file, overlap){
print(pheno.file)

  if(overlap == "overlap"){
  system(paste0("./mtg/mtg2 -fam ", paste0(geno.file, ".fam")," -bg ",
                paste0(geno.file,".grm.bin")," -pheno ", pheno.file,
                " -thread-num 10 -mod 2 -cove 1 -out ", output.file ))
  }

  if(overlap == "nooverlap"){
    system(paste0("./mtg/mtg2 -fam ", paste0(geno.file, ".fam")," -bg ",
                  paste0(geno.file,".grm.bin")," -pheno ", pheno.file,
                  " -thread-num 10 -mod 2 -out ", output.file ))
  }

}


CCSimuGCTA = function(geno.file,
                      pheno.file,
                      output.dir,
                      prev.Y = NULL,
                      prev.Z = NULL,
                      ncores = 1){

  pheno.file = paste0(output.dir, "outcome.pheno")

  pheno.outcome = data.table::fread(pheno.file)

  pheno.outcome %>%
    dplyr::select(FID,  IID) %>%
    fwrite( file = paste0(output.dir,"keep.txt"),
            quote = F, sep = "\t", col.names=T, na = "NA")

  grm = paste0(output.dir, "geno")

  system(paste0("./gcta64 --bfile ", geno.file, " --keep ", paste0(output.dir, "keep.txt"), " --make-grm-bin --out ", grm ) )

  output.file = paste0(output.dir, "GCTA.result")

  if(is.null(prev.Y)){

    system(paste0("./gcta64 --reml-bivar --grm ",
                  grm," --pheno ", pheno.file,
                  " --thread-num 20 --out ", output.file ))
  }else{

  system(paste0("./gcta64 --reml-bivar 1 2 --grm ",grm," --pheno ", pheno.file, " --reml-bivar-prevalence ",prev.Y, " ", prev.Z,
                " --thread-num 10 --out ", output.file ))
}
  file.remove(paste0(grm, ".grm.N.bin"))
  file.remove(paste0(grm, ".grm.bin"))

}
