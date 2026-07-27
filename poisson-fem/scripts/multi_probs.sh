#!/bin/bash

MODE=$1   # CPU or GPU
NP=$2     # total MPI ranks

if [[ -z "$MODE" || -z "$NP" ]]; then
  echo "Usage: ./multi_probs.sh CPU|GPU <num_processes>"
  exit 1
fi

#ENV SETUP
source /home/apps/spack/share/spack/setup-env.sh
spack load petsc@3.22.5 +cuda

export PETSC_DIR=$(spack location -i petsc@3.22.5)
export PETSC_ARCH=""
export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash

#COMPILE
rm -f multi_probs
mpif90 multi_probs.F90 -I$PETSC_DIR/include -L$PETSC_DIR/lib -lpetsc -o multi_probs

if [[ $? -ne 0 ]]; then
  echo "Compilation failed"
  exit 1
fi

#PETSc OPTIONS

if [[ "$MODE" == "CPU" ]]; then

  PETSC_OPTS=""

  echo "Running in CPU mode"
  echo "Total MPI ranks: $NP"

elif [[ "$MODE" == "GPU" ]]; then

  PETSC_OPTS="-mat_type aijcusparse -vec_type mpicuda -use_gpu_aware_mpi 1"

  echo "Running in GPU mode"
  echo "Total MPI ranks: $NP"

else
  echo "Invalid mode: use CPU or GPU"
  exit 1
fi

#RUN

mpirun -np $NP --bind-to none bash -c "
   export CUDA_VISIBLE_DEVICES=\$OMPI_COMM_WORLD_LOCAL_RANK
  ./multi_probs $PETSC_OPTS
"

