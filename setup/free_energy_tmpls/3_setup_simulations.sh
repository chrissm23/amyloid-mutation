#!/bin/bash
#
# Setup for simulations. Create folders and links to input files with the
# necessary modifications
#

source ../../../windows_res
states=($(printf "%s," "${windows[@]}"))
nstates=${#windows[@]}

basedir=../
top=$(pwd)
setup_dir=$(cd "$basedir"; pwd)

if [ \! -d simulations ]; then
  mkdir simulations
fi

cd simulations

for i in "${!windows[@]}"; do
  if [ \! -d ${windows[$i]} ]; then
    mkdir ${windows[$i]}
  fi

  sed -e "s/%L%/${windows[$i]}/;\
   s/%ne%/$nstates/;s/%states%/$states/" \
   $top/prod.tmpl > ${windows[$i]}/ti.in

  (
    cd ${windows[$i]}
    ln -sf $setup_dir/merged_topology.parm7 ti.parm7
    ln -sf $top/equilibration/${windows[$i]}/heat.rst7 ti.rst7
  )
done

cd $top
