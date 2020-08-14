#phase with SHAPEIT - pop and read back.
########################################

#split by scafs
while read scafs
do
vcftools --vcf pseudo_qual_dust_indel_fixedmissing_nodups.vcf --recode --recode-INFO-all --max-missing .9 --out ${scafs}_missing10p --chr ${scafs} & 
done < scaf_list

mv /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/filtering/Scaffold_*_missing10p.recode.vcf /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/

#copy map files
cp /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/ldhat/LDhat_workflow/amaranth/*_finalformat_recombrate.txt /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing

#get scaf list
cp /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/filtering/scaf_list /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing

#get SNP positions
while read scafs
do
grep -v "#" ${scafs}_missing10p.recode.vcf | awk '{print $2}' > ${scafs}.pos
done < scaf_list

#check for duplicates in map file
while read scafs
do
uniq ${scafs}_finalformat_recombrate.txt > ${scafs}_finalformat_recombrate_uniq.txt
done < scaf_list

#impute recombination rates between SNPs in our VCF
#arg1 = map file, arg2 = listofpos, arg3 = output name
while read scafs
do
Rscript imputemap.R ${scafs}_finalformat_recombrate_uniq.txt ${scafs}.pos ${scafs}_shapeit.map 
done < scaf_list

#remove unmatched position in vcf
while read scafs
do
#awk -v var=${scafs} '{print var "\t" $1 }' ${scafs}_shapeit.map > ${scafs}_map.pos
awk '{if ($2 <= 0) print $1 " " 0.001 " " $3; else print $1 " " $2 " " $3}' ${scafs}_shapeit.map | sed '1 i\pposition rrate gposition' > ${scafs}_shapeit2.map
#vcftools --vcf ${scafs}_missing10p.recode.vcf --recode --recode-INFO-all --positions ${scafs}_map.pos --out ${scafs}_missing_mapmatched &
done < scaf_list


#run shape it
while read scafs
do
shapeit --input-vcf ${scafs}_missing_mapmatched.recode.vcf --input-map ${scafs}_shapeit2.map --window 0.5 --effective-size 500000 --thread 10 --output-max ${scaf}_shapeit  &
done < scaf_list

#############################################################
