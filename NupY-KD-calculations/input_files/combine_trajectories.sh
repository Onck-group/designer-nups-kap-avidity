#!/bin/bash

#assumes 20 sets of replicas
gmx trjcat -f replica_*/traj_comp.xtc -o traj_combined.xtc -settime << EOF
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
c
EOF
