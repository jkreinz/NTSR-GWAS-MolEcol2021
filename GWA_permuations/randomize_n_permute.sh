##########################################
#run this script using parallel_permute.sh
##########################################
#within population GWA permutations (geographic region population embedded but commented out)

#set up variables
rep=$1
mainpath='/ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/filtering/gwas/'

#remove files if they already exist since appending loops
rm -R randpheno_bypop_${rep}
mkdir randpheno_bypop_${rep}
cd randpheno_bypop_${rep}

#copy over intermediate files to run in each directory
paste -d" " $mainpath/bypop_forpheno.fam $mainpath/covariates > intermed

#permute by population
while read pop
do
sed 's/N./N/g' $mainpath/intermed | awk -v pop=${pop} '$1 == pop' | shuf | cut -d$'\t' -f6,7,8 >> randpheno_bypop_${rep} #jointly randomizes resistance and covariates by population, Natural (N) samples treated as one population
done < $mainpath/pops

#permute by geographic region
#rm randpheno_byregion_${rep}
#paste <(grep "^[4-9]" $mainpath/intermed | cut -d$'\t' -f-2) <( grep "^[4-9]" $mainpath/intermed | shuf | cut -d$'\t' -f6,7,8) >> randpheno_byregion_${rep}
#paste <(grep "^[0-9][0-9]" $mainpath/intermed | cut -d$'\t' -f-2) <( grep "^[0-9][0-9]" $mainpath/intermed | shuf | cut -d$'\t' -f6,7,8) >> randpheno_byregion_${rep}
#paste <( grep -E "^B|^D|^E|^F|^G|^H|^J" $mainpath/intermed | cut -d$'\t' -f-2) <(grep -E "^B|^D|^E|^F|^G|^H|^J" $mainpath/intermed | shuf | cut -d$'\t' -f6,7,8) >> randpheno_byregion_${rep}
#paste <(grep -E "^K" $mainpath/intermed | cut -d$'\t' -f-2) <(grep -E "^K" $mainpath/intermed | shuf | cut -d$'\t' -f6,7,8 ) >> randpheno_byregion_${rep}
#paste <(grep -E "^N[0-9]" $mainpath/intermed |  cut -d$'\t' -f-2) <(grep  -E "^N[0-9]" $mainpath/intermed | shuf | cut -d$'\t' -f6,7,8) >> randpheno_byregion_${rep}

#merge randomized phenos with genotype labels
##fam file for pops
awk '{print NR " " $0}' intermed | cut -d$' ' -f1,2,3 > reorder #put back to original order as in .fam
paste -d" " <(while read pop; do cat intermed | awk -v pop=${pop} '$1 == pop' | cut -d$' ' -f1-2; done < $mainpath/../pops_22) randpheno_bypop_${rep} > randpheno_bypop_${rep}_wsampnames
paste -d" " <(sort -k2,2 -k3,3 reorder) <(sort -k1,1 -k2,2 randpheno_bypop_${rep}_wsampnames | sed 's/^ //g') | sort -nk1,1 | cut -d$' ' -f4- | awk '{print $1 "\t" $2 "\t" 0 "\t" 0 "\t" 0 "\t" $3 "\t" $4 "\t" $5}' > forgwas_nodups_tryagain_updated_indsexcluded_nonbin_perm${rep}.fam
cut -d$'\t' -f7- forgwas_nodups_tryagain_updated_indsexcluded_nonbin_perm${rep}.fam > covs

##fam file for geographic region permutations
#awk '{print NR "\t" $0}' intermed | cut -d$'\t' -f1,2,3 > reorder #put back to original order as in .fam
#paste <(sort -k2,2 -k3,3 reorder) <(sort -k1,1 -k2,2 randpheno_byregion_${rep}) | sort -nk1,1 | cut -d$'\t' -f4- | awk '{print $1 "\t" $2 "\t" 0 "\t" 0 "\t" 0 "\t" $3 "\t" $4 "\t" $5}' > forgwas_nodups_tryagain_updated_indsexcluded_perm${rep}.fam

#copy over GEMMA input files into working directory & name fam to match bam/bim files
cp forgwas_nodups_tryagain_updated_indsexcluded_perm${rep}.fam forgwas_nodups_tryagain_updated_indsexcluded.fam
cp $mainpath/forgwas_nodups_tryagain_updated_indsexcluded.b?? ./

#run permuted GEMMA GWA with covariates
/ohta/julia.kreiner/software/gemma-0.98.1-linux-static -bfile forgwas_nodups_tryagain_updated_indsexcluded -lmm 4 -k $mainpath/../output/forgwas_nodups_tryagain_updated_indsexcluded.cXX.txt -o ../forgwas_nodups_tryagain_updated_indsexcluded_popperm_covariates_${rep} -c covs
rm forgwas_nodups_tryagain_updated_indsexcluded.b*

cd ../
