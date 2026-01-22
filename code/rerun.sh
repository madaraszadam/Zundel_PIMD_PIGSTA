export PATH=/home/madarasz/simulations/CP2K/Zundel_cation/PIMD/PILE/T1.67K/code:$PATH

export CP2K_DATA_DIR="/home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1/data"

export OMP_NUM_THREADS=1

#setenv CP2K_DATA_DIR "/home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1/data"

bname=Zundel_cation_PIMD

filename=filt_$bname

ride=400
previ=400026
natoms=$(head -n 1 filt_Zundel_cation_PIMD-pos-1-1.xyz )
#nframe=40
Pbead=NUMBER_OF_BEADS

first=$(($(head $filename-pos-1-1.xyz | sed -n 2,2p | sed 's@^[^0-9]*\([0-9]\+\).*@\1@')))

echo $first

last=$(($(tail -n $(($natoms+2)) $filename-pos-1-1.xyz | sed -n 2,2p | sed 's@^[^0-9]*\([0-9]\+\).*@\1@')))

echo $last

#rm -f samp_$filename-*-1.xyz

#rm -f samp_$bname-*-1.xyz

for i in `seq 1 $Pbead`; do

  for j in $(seq $first $ride $last) ; do grep " $j," -m 1 -A $natoms -B 1 $filename-pos-$i-1.xyz >> samp_$filename-pos-$i-1.xyz; done
  for j in $(seq $first $ride $last) ; do grep " $j," -m 1 -A $natoms -B 1 $filename-vel-$i-1.xyz >> samp_$filename-vel-$i-1.xyz; done
  for j in $(seq $first $ride $last) ; do grep " $j," -m 1 -A $natoms -B 1 $filename-totforce-$i-1.xyz >> samp_$filename-totforce-$i-1.xyz; done

  for j in $(seq $first $ride $last) ; do grep " $j," -m 1 -A $natoms -B 1 ../prodrun1/$bname-pos-$i-1.xyz >> samp_$bname-pos-$i-1.xyz; done
  for j in $(seq $first $ride $last) ; do grep " $j," -m 1 -A $natoms -B 1 ../prodrun1/$bname-totforce-$i-1.xyz >> samp_$bname-totforce-$i-1.xyz; done

#  split -l $((($natoms+2))) -d -a 6 samp_$filename-pos-$i-1.xyz fp$i"_"
#  split -l $((($natoms+2))) -d -a 6 samp_$filename-vel-$i-1.xyz fv$i"_"

done

nframe=$(grep -c "E =" samp_$filename-pos-1-1.xyz)
echo $nframe

#sep='\"\ \"'
#echo $sep
#formp=\'\{print\ $(eval echo 'bohr\*\$'{2..$(($Pbead*4-2))..4}$sep 'bohr\*\$'{3..$(($Pbead*4-1))..4}$sep 'bohr\*\$'{4..$(($Pbead*4))..4}$sep)\}\'
#formv=\'\{print\ $(eval echo '\$'{2..$(($Pbead*4-2))..4}$sep '\$'{3..$(($Pbead*4-1))..4}$sep '\$'{4..$(($Pbead*4))..4}$sep)\}\'
#echo "$formp"


listfr=$(eval echo {000000..$(($nframe-1))})


head -n 1 ../prodrun1/Zundel_cation_PIMD-energy-1.dat > energy.dat

rm recalc-force-*-1.xyz

iter=$first

for j in $listfr; do

#  list=$(for j in `seq 1 $Pbead`; do printf "fp$j_$(($i-1)) " ; done)

#  listp=$(eval echo fp{1..$Pbead}_$j)
#  listv=$(eval echo fv{1..$Pbead}_$j)
#  sed -i -e 1,2d $listp $listv

#  echo $listp
#  paste $listp | awk -v bohr=1.8897259886 '{ print bohr*$2" "bohr*$6" "bohr*$10" "bohr*$14" "bohr*$3" "bohr*$7" "bohr*$11" "bohr*$15" "bohr*$4" "bohr*$8" "bohr*$12" "bohr*$16"\\"}' > bead_pos.xyz
#  paste $listp | tail -n $natoms | eval awk -v bohr=1.8897259886 "$formp" > bead_pos.xyz
#  paste $listp | tail -n $natoms | eval awk "$formv" > bead_pos.xyz

#  sed -i 's/ /\n/g' bead_pos.xyz
#  sed -i '/^$/d' bead_pos.xyz
#  awk -v factor=1.8897259886 '{ printf "%.10f\n", $1*factor}' bead_pos.xyz > bead_pos_bohr.xyz
#  mv bead_pos_bohr.xyz bead_pos.xyz
#  sed -i -e 's/$/ \\/' -e '$s/\\$//' bead_pos.xyz



#  echo $listv

#  paste $listv | tail -n $natoms | eval awk "$formv" > bead_vel.xyz

#  sed -i 's/ /\n/g' bead_vel.xyz
#  sed -i '/^$/d' bead_vel.xyz
#  sed -i -e 's/$/ \\/' -e '$s/\\$//' bead_vel.xyz

#  truncate -s -2 bead_vel.xyz

#  rm recalc_filt*

  echo "samp_$filename-pos-" > input.txt
  echo $Pbead >> input.txt
  xyz-cp2k  < input.txt

  mv cp2k_samp_$filename-pos-1-1.xyz bead_pos.xyz

  echo "samp_$filename-vel-" > input.txt
  echo $Pbead >> input.txt
  xyz-cp2k  < input.txt

  mv cp2k_samp_$filename-vel-1-1.xyz bead_vel.xyz

  rm recalc_filt-energy-1.dat recalc_filt-r-0.out

  rm recalc_filt-force-*-1.xyz

#  cp ../../recalc_filt.inp .

  sed -i "s/ITERATION  $previ/ITERATION  $iter/g" recalc_filt.inp

  echo "$iter"

  /home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1-Linux-gnu-x86_64.ssmp recalc_filt.inp

  tail -n 1 recalc_filt-energy-1.dat >> energy.dat

  for i in `seq 1 $Pbead`; do

    cat recalc_filt-force-$i-1.xyz >> recalc-force-$i-1.xyz

  done

  previ=$iter

  iter=$(($iter+$ride))

  rm bead_pos.xyz

  rm bead_vel.xyz

  rm $listp

  rm $listv

done

rm fp*_*
rm fv*_*

