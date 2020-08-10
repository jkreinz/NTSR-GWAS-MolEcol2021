i=$1

rm -R Scaf${i}_IHS
mkdir Scaf${i}_IHS

#while read hits
#do
selscan --hap Scaffold_${i}_selscan.haps --map  Scaffold_${i}_selscan.map --ihs --ihs-detail --out Scaf${i}_IHS/Scaf${i} --maf 0.01
#done < Scaf${i}_hits
