#!/bin/bash
#
# Run simulations to correct for the density of the solvation box at lambda=0
#
#SBATCH --job-name=%r%_%rn%_solv_%pairID%
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber/20.0/20.4

prmtop=../merged_topology.parm7
inpcrd=../merged_topology.rst7

echo "Minimizing..."
mpirun -n 8 pmemd.MPI -i min.in -p $prmtop -c $inpcrd -ref $inpcrd \
  -O -o min.out -e min.en -inf min.info -r min.rst7

echo "Heating..."
mpirun -n 8 pmemd.MPI -i heat.in -p $prmtop -c min.rst7 -ref $inpcrd \
  -O -o heat.out -e press.en -inf heat.info -r heat.rst7 -x heat.nc

echo "Pressuring..."
mpirun -n 8 pmemd.MPI -i press.in -p $prmtop -c heat.rst7 -ref heat.rst7 \
  -O -o press.out -e press.en -inf press.info -r press.rst7 -x press.nc

echo "Done"
