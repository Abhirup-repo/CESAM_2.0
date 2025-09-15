#!/bin/bash -ex

module purge
module list
module load git
#module load openmpi/4.1.4-static-intel22
#module load intel/23.2.1
module load openmpi/5.0.5-static-intel24
module load intel/24.0.1
module load cdo/2.2.2-gcc-12.2.0
