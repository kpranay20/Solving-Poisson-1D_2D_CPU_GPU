#!/bin/bash

#Terminal excecution of this jobscriopt:
#   ./petsc_gpu.sh CPU 64 : If want to run on CPU with 64 mpi ranks 
#   ./petsc_gpu.sh GPU 1 : If want to run on GPU with 1 node

MODE=$1
NP=$2

if [[ -z "$MODE" || -z "$NP" ]]; then
  echo "Usage: ./run.sh CPU|GPU <num_processes>"
  exit 1
fi

#Environment loading
source /home/apps/spack/share/spack/setup-env.sh
spack load petsc@3.22.5 +cuda

#PETSc variables declaration(uses Spack prefix automatically)
export PETSC_DIR=$(spack location -i petsc@3.22.5)
export PETSC_ARCH=""

export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash

#Cleaning and compiling
rm -f poisson

mpif90 petsc_gpu.F90 -I$PETSC_DIR/include -L$PETSC_DIR/lib -lpetsc -o poisson

#---------------------------------------------------------------------------------

if [[ $? -ne 0 ]]; then
  echo "Compilation failed"
  exit 1
fi

# PETSc runtime options
if [[ "$MODE" == "CPU" ]]; then
  PETSC_OPTS=""
elif [[ "$MODE" == "GPU" ]]; then
  PETSC_OPTS="
    -mat_type aijcusparse
    -vec_type mpicuda
    -use_gpu_aware_mpi 1
    -log_view
  "
else
  echo "Invalid mode: use CPU or GPU"
  exit 1
fi

#Running based on mode (CPU/GPU) 
echo "Mode       : $MODE"
echo "MPI ranks  : $NP"

mpirun -np $NP ./poisson $PETSC_OPTS

