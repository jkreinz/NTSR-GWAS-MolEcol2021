seq 1 250 | parallel -j 20 "bash randomize_n_permute.sh {}" #first produced randomized sample files & GWAS outputs

#then used randomized files to calculate Fst relative to randomized res-sus comparisons for that iteration
for i in {1..250}
do
	awk '$6 >= 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > res_inds_${i} #list interations resistant inds
	awk '$6 < 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > sus_inds_${i} #list of interations sus inds

	vcftools --vcf TSR_permutedsigsites.recode.vcf --weir-fst-pop TSR_rep_${i}/res_inds --weir-fst-pop res_inds_${i} --weir-fst-pop sus_inds_${i} --out TSR_permutedsigsites_RvS_${i}
	#calculate weir cockerham fst
done

cat *fst | grep -v "CHROM" > TSR_permutedsigsites_RvS.weir.fst #combine fst estimates for all sites across all iterations to estimate empirical null expectation
