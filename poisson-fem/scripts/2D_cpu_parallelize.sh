#Depending on the PETSc installation, there might be ton of errors, make sure you know the PETSC_DIR and PETSC_ARCH locations. All the errors revolve around these two.

#In our case, the 2-GPU node machine, a PETSc env was created so this worked:
source /home/apps/spack/share/spack/setup-env.sh
spack load petsc@3.22.5 +cuda
export PETSC_DIR=$(spack location -i petsc@3.22.5)
export PETSC_ARCH=""
export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash

#Common for all:
mpif90 2D_cpu_parallelize.F90 -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include -L$PETSC_DIR/$PETSC_ARCH/lib -lpetsc -o petsc_parallel
mpirun -np 4 ./petsc_parallel   #This will run the code on 4 CPU nodes.•
