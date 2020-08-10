i=$1

rm -R Scaf${i}_IHS
mkdir Scaf${i}_IHS

while read hits
do
selscan --hap Scaffold_${i}_selscan.haps --map  Scaffold_${i}_selscan.map --ehh ${hits} --out Scaf${i}_EHH/Scaf${i} --maf 0.01
done < Scaf${i}_hits
