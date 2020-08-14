##########################################
#run this script using parallel_permute.sh
##########################################


rep=$1

#organize directories
rm -R TSR_rep_${rep}
mkdir TSR_rep_${rep}
cd TSR_rep_${rep}

#make files for GWAS iteration
paste -d" " ../../forgwas_nodups_tryagain_updated_indsexcluded.fam ../../covariates | shuf | cut -d" " -f6,7,8 > randpheno_${rep} #keep resistance level, copy number, TSR status
paste -d" " ../ind_order <(cut -d" " -f1 randpheno_${rep}) > forgwas_nodups_tryagain_updated_indsexcluded_perm${rep}.fam #ind order

#keep track
cp forgwas_nodups_tryagain_updated_indsexcluded_perm${rep}.fam forgwas_nodups_tryagain_updated_indsexcluded.fam
cut -d" " -f2 randpheno_${rep} > TSR_${rep}
cp ../forgwas_nodups_tryagain_updated_indsexcluded.b* ./

#run GEMMA on randomized files
/ohta/julia.kreiner/software/gemma-0.98.1-linux-static -bfile forgwas_nodups_tryagain_updated_indsexcluded -lmm 4 -k ../../output/forgwas_nodups_tryagain_updated_indsexcluded.cXX.txt -o ../forgwas_nodups_tryagain_updated_indsexcluded_TSR_${rep} -c TSR_${rep}

#rm forgwas_nodups_tryagain_updated_indsexcluded.b*

cd ../
