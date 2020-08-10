seq 1 250 | parallel -j 20 "bash randomize_n_permute.sh {}"
