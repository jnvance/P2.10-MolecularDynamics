#!/bin/bash
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=20
# PBS -l nodes=3:ppn=20
#PBS -q reserved3
#PBS -T flush_cache

## For Ulysses
module load testing; 
module load openmpi/2.0.0/gnu/6.2.0

cd $PBS_O_WORKDIR
make clean
make

results=$PBS_O_WORKDIR/data_$PBS_JOBID
mkdir $results

timings=$results/timings_$PBS_JOBID.txt
out=$results/out_$PBS_JOBID.txt


echo "Timings: $timings"
echo "Lattice: $lattice"

cd $results;

# ../../input/gen_input_noexchange.py
../../input/gen_input.py

lattice=15
../../input/lattice $lattice > crystal.xyz 

for nprocs in $(seq 16 16)
do
    telapsed=$(( (/usr/bin/time -p -f%e mpirun -np $nprocs ../simplemd.x *.params) 1>>$out)  2>&1)
    echo $lattice $nprocs $telapsed >> $timings
    $out
done 

lattice=18
../../input/lattice $lattice > crystal.xyz 

for nprocs in $(seq 16 16)
do
    telapsed=$(( (/usr/bin/time -p -f%e mpirun -np $nprocs ../simplemd.x *.params) 1>>$out)  2>&1)
    echo $lattice $nprocs $telapsed >> $timings
    $out
done 

