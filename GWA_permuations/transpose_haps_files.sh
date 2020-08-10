i=$1

datamash --no-strict transpose <  <(cut -d" " -f6- ~/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_${i}_shapeit.haps | tr " " "\t") > Scaffold_${i}_selscan.haps
