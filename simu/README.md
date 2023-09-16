This folder contains the code for the simulation: 

Here is the example for the GWAS simulation 



Simulate the phenotype 

```{bash}
cd /home/wangjq/Linear-based/GWAS_Simu_2023/
filename="/home/wangjq/Linear-based/GWAS_Simu_2023/pheno/Pheno_CON.R"
genofile="/home/wangjq/PMBB/prepare/PMBB_geno_EU"
coef='normal'
overlap='overlap'
n_sample=20000

args_list=("0.001 0.001 0.2 0.001" "0.002 0.002 0.2 0.002" "0.01 0.01 0.2 0.01" "0.1 0.0025 0.2 0.0025" "1 0.0005 0.2 0.0005" "0.1 0.1 0.2 0.0025" "0.05 0.05 0.2 0.05")

for coef in 'LDAK' 'LDAKs3' 'INTER' #'normal' 
do
for args in "${args_list[@]}"; do
    bsub -n 5 -M 60240 -R "rusage [mem=60240]" -o ./temp/${current_time}_pheno_output.txt "Rscript ${filename} ${genofile}  ${args} ${coef} ${overlap} ${n_sample}"
done
done

```


## Main function 

```{bash}
cd /home/wangjq/Linear-based/GWAS_Simu_2023/
module load R
filename='./main/Main2.R'
coef='normal'
overlap='overlap'
args_list=("0.001 0.001 0.2 0.001" "0.002 0.002 0.2 0.002" "0.01 0.01 0.2 0.01" "0.1 0.0025 0.2 0.0025" "1 0.0005 0.2 0.0005" "0.1 0.1 0.2 0.0025" "0.05 0.05 0.2 0.05")

args_list=("0.001 0.001 0.2 0.001" "0.01 0.01 0.2 0.01")


for args in "${args_list[@]}"; do
  bsub -n 1 -J "jobname[1-10]%10" -M 60240 -R "rusage [mem=60240]" -o test.out "Rscript ${filename} --args ${args} ${coef} ${overlap}"
done

filename='./MainCode/MainCompare.R'


```
