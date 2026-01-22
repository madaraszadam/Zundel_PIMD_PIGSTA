cat ./SEED*/Zundel_cation_PIMD-energy.dat > PIMD_energies.txt
tail -q -n +2 ./SEED*/energy.dat > filtered_energies.txt
paste PIMD_energies.txt filtered_energies.txt > all_energies.txt
