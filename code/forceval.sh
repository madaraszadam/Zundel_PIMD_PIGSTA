grep "E =" -A 7 --no-group-separator ./SEED*/recalc-totforce.xyz | grep -v "E =" > onlyrecalcforces.txt
grep "E =" -A 7 --no-group-separator ./SEED*/samp_filt_Zundel_cation_PIMD-totforce.xyz | grep -v "E =" > only_filt_totforces.txt
paste onlyrecalcforces.txt only_filt_totforces.txt > comp_forces.txt
