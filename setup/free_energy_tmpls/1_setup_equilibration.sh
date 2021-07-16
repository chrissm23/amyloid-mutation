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

if [ \! -d equilibration ]; then
  mkdir equilibration
fi

cd equilibration

for i in "${!windows[@]}"; do
  if [ \! -d ${windows[$i]} ]; then
    mkdir ${windows[$i]}
  fi

  sed -e "s/%L%/${windows[$i]}/" \
   $top/heat.tmpl > ${windows[$i]}/heat.in

  (
    cd ${windows[$i]}
    ln -sf $setup_dir/merged_topology.parm7 ti.parm7
    if [ $i == 0 ]; then
      ln -sf $setup_dir/minimization/min.rst7 ti.rst7
    fi
  )
done

cd $top
