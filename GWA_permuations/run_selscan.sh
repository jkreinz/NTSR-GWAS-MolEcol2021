seq 1 16 | parallel -j 16 "bash transpose_haps_files.sh {}"

#reformat map files for selscan
for i in {1..16}
do
awk -v var=$i '{print var "\t" "Scaf" var "_" $1 "\t" $3 "\t" $1}'  /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_${i}_shapeit2.map | tail -n +2 > Scaffold_${i}_selscan.map
done

#get seperate files of sig GWA hits by scaffold
cut -d"_" -f1 ../TSRsighits | uniq > scafs
while read scafs; do grep $scafs ../TSRsighits > ${scafs}_hits; done < scafs &

#run selscan in parallel across scaffolds
seq 1 16 | parallel -j 16 "bash get_EHH.sh {}"
seq 1 16 | parallel -j 16 "bash get_IHS.sh {}"
#and then do the same for permuted SNPs


