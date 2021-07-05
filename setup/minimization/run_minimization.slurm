#!/bin/bash
#
# Run minimization at lambda=0 for later equilibration
#
#SBATCH --job-name=%r%_%rn%_minimization_%pairID%
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber

prmtop=../merged_topology.parm7
inpcrd=../solvation_correction/press.rst7

echo "Minimizing..."
mpirun -n 8 pmemd.MPI -i minimization.in -p $prmtop -c $inpcrd -ref $inpcrd\
  -O -o min.out -e min.en -inf min.info -r min.rst7

echo "Done"
