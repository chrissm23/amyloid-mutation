#!/bin/bash
#
# Submit simulations to cow. 
# REMEMBER TO ADD WINDOW SET TO JOB NAME AND CHANGE WINDOWS
# REMEMBER TO ADD -ref FILE WHEN USING RESTRAINTS
#
#SBATCH --job-name=%r%_%rn%_prod_%dir%_
#SBATCH --time=4-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber/20.0/20.4

windows=(0.0 0.1 0.2 0.3 0.4 0.5)

cd simulations

for i in "${!windows[@]}"; do
  cd ${windows[$i]}
  echo "lambda = ${windows[$i]}"

  echo "Running MD..."
  pmemd.cuda -i ti.in -c ti.rst7 -p ti.parm7 \
  -O -o ti001.out -inf ti001.info -e ti001.en -r ti001.rst7 -x ti001.nc \
  -l ti001.log

  cd ..
done

cd ..

echo "Done"
