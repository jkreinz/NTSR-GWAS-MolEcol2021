seq 1 250 | parallel -j 20 "bash randomize_n_permute.sh {}"

for i in {1..250}
do
	awk '$6 >= 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > res_inds_${i}
	awk '$6 < 2' TSR_rep_${i}/forgwas_nodups_tryagain_updated_indsexcluded.fam | cut -d " " -f 1 > sus_inds_${i}

	vcftools --vcf TSR_permutedsigsites.recode.vcf --weir-fst-pop TSR_rep_${i}/res_inds --weir-fst-pop res_inds_${i} --weir-fst-pop sus_inds_${i} --out TSR_permutedsigsites_RvS_${i}
done

cat *fst | grep -v "CHROM" > TSR_permutedsigsites_RvS.weir.fst
