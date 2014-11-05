#!/bin/bash
#PBS -S /bin/bash
#PBS -A fgkstar
#PBS -l nodes=6:ppn=4
#PBS -l walltime=00:01:00
#PBS -m bea -M aivaz@dps.uminho.pt
#PBS -j oe
cd $PBS_O_WORKDIR
module load comp/openmpi/gcc4
mpipbsexec -np 24 -machinefile $PBS_NODEFILE ./pswarm > output.txt 2>&1
