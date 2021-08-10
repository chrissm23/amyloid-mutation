#!/bin/bash
#
# Combine trajectories from all lambada windows, autoimage and strip water
#

cpptraj -p ../merged_topology.parm7 <<_EOF

trajin simulations/*/ti001.nc 1 last 10
autoimage
strip :WAT
trajout full_traj_nowater.nc
go
exit

_EOF

cpptraj -p ../merged_topology.parm7 <<_EOF

parmstrip :WAT
parmwrite out merged_topology_nowater.parm7
go
exit

_EOF
