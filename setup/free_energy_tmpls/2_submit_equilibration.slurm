#!/bin/bash
#
# Submit simulations to cow. REMEMBER TO ADD JOB NAME
#
#SBATCH --job-name=%r%_%rn%_equil_%dir%
#SBATCH --time=4-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber/20.0/20.4

source ../../../windows_res

cd equilibration

for i in "${!windows[@]}"; do
  cd ${windows[$i]}
  echo "lambda = ${windows[$i]}"

  if [ $i == 0 ]; then
    inpcrd=./ti.rst7
  else
    inpcrd=../${windows[$i-1]}/heat.rst7
  fi
  echo "Heating..."
  pmemd.cuda -i heat.in -c $inpcrd -ref $inpcrd -p ti.parm7 \
    -O -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log

  cd ..
done

cd ..

echo "Done"
