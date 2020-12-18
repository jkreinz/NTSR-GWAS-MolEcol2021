seq 1 250 | parallel -j 20 "bash randomize_n_permute.sh {}" #first produced randomized sample files & GWAS outputs

#then used randomized files to calculate Fst relative to randomized res-sus comparisons for that iteration
for i in {1..250}
do
	awk '$6 >= 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > res_inds_${i} #list interations resistant inds
	awk '$6 < 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > sus_inds_${i} #list of interations sus inds

	vcftools --vcf TSR_permutedsigsites.recode.vcf --weir-fst-pop res_inds_${i} --weir-fst-pop sus_inds_${i} --out TSR_permutedsigsites_RvS_${i}
	#calculate weir cockerham fst
done

cat *fst | grep -v "CHROM" > TSR_permutedsigsites_RvS.weir.fst #combine fst estimates for all sites across all iterations to estimate empirical null expectation




for i in {1..250}; do awk '$6 > 0' ../Binarypop_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | awk '{print $1 "\t" $2}' | sed -e's/K\t/K/g' | sed -e's/D\t/D/g' | sed -e's/B\t/B/g' | sed -e's/H\t/H/g' | sed -e's/E\t/E/g' | sed -e's/F\t/F/g' | sed -e's/G\t/G/g' | sed -e's/J\t/J/g' | tr "\t" "_" > res_inds_${i}; awk '$6 < 1' ../Binarypop_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | awk '{print $1 "\t" $2}' | sed -e's/K\t/K/g' | sed -e's/D\t/D/g' | sed -e's/B\t/B/g' | sed -e's/H\t/H/g' | sed -e's/E\t/E/g' | sed -e's/F\t/F/g' | sed -e's/G\t/G/g' | sed -e's/J\t/J/g' | tr "\t" "_" > sus_inds_${i}; grep "^$i    " ../binary_FDR05_sigsnps.txt | awk '{print $2 "\t" $4}' > rep_${i}_pos; vcftools --vcf ../binary_FDR05_sigsnps.recode.vcf --weir-fst-pop res_inds_${i} --weir-fst-pop sus_inds_${i} --out binary_FDR05_permutedsigsites_RvS_${i} --positions rep_${i}_pos; done
