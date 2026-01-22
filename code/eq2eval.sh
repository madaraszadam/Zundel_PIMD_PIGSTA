Pbead=4
filename=Zundel_cation_PIMD

export LD_LIBRARY_PATH=/home/rokoba/amber/lib:$LD_LIBRARY_PATH
export CP2K_DATA_DIR="/home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1/data"
export OMP_NUM_THREADS=1

ulimit -n 4096

tmpdir=$(pwd)
export PATH=$tmpdir/../code:$PATH
mkdir -p /state/partition1/madarasz$tmpdir
cd /state/partition1/madarasz$tmpdir

echo $1 
rm -r SEED$1
mkdir SEED$1
cd SEED$1
mkdir equilibration
cd equilibration
cp $tmpdir/$filename.inp .
cp $tmpdir/opt_geom.xyz .
sed -i "s/THERMOSTAT_SEED 1/THERMOSTAT_SEED $1/g" $filename.inp 
sed -i "s/NUMBER_OF_BEADS/$Pbead/g" $filename.inp
/home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1-Linux-gnu-x86_64.ssmp $filename.inp
cd ..

mkdir prodrun1
cd prodrun1

cp $tmpdir/../kernel_library/P"$Pbead"_T1.67K_dt0.25fs/kernel_function.dat kernel.dat
cp ../equilibration/$filename-1.restart .
cp ../equilibration/opt_geom.xyz .

klen=$(cat kernel.dat | wc -l)
slen=$((2*$klen-2))


sed -i "s/VELOCITIES  OFF/VELOCITIES  ON/g" $filename-1.restart
sed -i "s/FORCES  OFF/FORCES  ON/g" $filename-1.restart
sed -i "s/FORMAT  ATOMIC/FORMAT  XYZ/g" $filename-1.restart
sed -i "s/PINT  400/PINT  1/g" $filename-1.restart
sed -i "s/NUM_STEPS  40000/NUM_STEPS  $slen/g" $filename-1.restart
sed -i "s/TAU     9.9999999999999989E+02/TAU     0/g" $filename-1.restart
sed -i "s/LAMBDA     5.0000000000000000E-01/LAMBDA     0.0/g" $filename-1.restart

/home/madarasz/program_dev/cp2k-2023.1/cp2k-2023.1-Linux-gnu-x86_64.ssmp  $filename-1.restart

echo "$filename-pos-" > input.txt
echo $Pbead >> input.txt
echo "1.67" >> input.txt

spring_force_calc < input.txt 

#for i in `seq 1 $Pbead`; do
for i in `seq 1 1`; do

  echo "$filename-force-$i-1.xyz" > input.txt
  echo "spf_$filename-pos-$i-1.xyz" >> input.txt
  echo "$filename-totforce-$i-1.xyz" >> input.txt

  sum_forces < input.txt 

done

cd ..

mkdir filter1
cd prodrun1
#  cp ../prodrun1/$filename*.xyz .

for i in `seq 1 $Pbead`; do

  echo "$filename-pos-$i-1.xyz" > input.txt
  echo "0.00025" >> input.txt
  echo "1.67" >> input.txt

  filteronly < input.txt 

  sed -i "s/pos/vel/g" input.txt
  filteronly < input.txt

  if [ $i -eq 1 ]
  then

    sed -i "s/vel/totforce/g" input.txt
    filteronly < input.txt

  fi

done

mv filt_*.xyz ../filter1
#  rm $filename*.xyz
cd ..

echo $1 
cd filter1
cp $tmpdir/../code/rerun.sh .
cp /$tmpdir/../code/recalc_filt.inp .
cp $tmpdir/opt_geom.xyz .
sed -i "s/NUMBER_OF_BEADS/$Pbead/g" rerun.sh  
sed -i "s/NUMBER_OF_BEADS/$Pbead/g" recalc_filt.inp
./rerun.sh > /dev/null

echo "samp_filt_$filename-pos-" > input.txt
echo $Pbead >> input.txt
echo "1.67" >> input.txt

spring_force_calc < input.txt

for i in `seq 1 1`; do

  echo "recalc-force-$i-1.xyz" > input.txt
  echo "spf_samp_filt_$filename-pos-$i-1.xyz" >> input.txt
  echo "recalc-totforce-$i-1.xyz" >> input.txt
  sum_forces < input.txt

done

rm *.txt

echo "samp_filt_$filename-pos-" > g-a.inp
echo $Pbead >> g-a.inp
bead_anal < g-a.inp > g-a.out

rename txt filt *.txt

echo "samp_$filename-pos-" > g-a.inp
echo $Pbead >> g-a.inp
bead_anal < g-a.inp > g-a.out

rename txt pimd *.txt

echo "cleaning up"

rename -- -1-1 "" *-1-1.xyz

mv recalc_filt-centroid-gyr-1.dat recalc-centroid-gyr.dat

sed "$(($klen))q;d" ../prodrun1/$filename-centroid-gyr-1.dat > $filename-centroid-gyr.dat

sed "$(($klen+1))q;d" ../prodrun1/$filename-energy-1.dat > $filename-energy.dat

mv ../equilibration/$filename-1.restart .

rm g-a.* recalc_f*.* rerun.sh filt_* *1.xyz

cd ..

rm -r equilibration

rm -r prodrun1

mv filter1 $tmpdir/SEED$1

cd ..

rm -r SEED$1

