# Finite Element Solution of the Poisson Equation  

## Overview

This repository documents the development of a **finite element solver for the Poisson equation**, starting from simple CPU-based implementations and extending to parallel and later GPU-accelerated solvers using PETSc.

The same mathematical problem is solved across multiple programming languages and hardware architectures in order to verify the numerical correctness and study performance scaling in MPI and GPU solvers.


## Problem Description

We consider the Poisson equation

$$\nabla^2 u(\mathbf{x}) = f(\mathbf{x})$$

on one- and two-dimensional domains with Dirichlet boundary conditions. The problem is discretized using the finite element method (FEM) with low-order basis functions and explicit matrix assembly. Full mathematical formulation and discretization details are provided in the documentation.



## Implementation Workflow

The solver is developed in the following manner:

1. **Python (NumPy, CPU)**

2. **Python (CuPy, GPU)**: GPU Implementation in Python (Without MPI)

3. **Fortran + LAPACK (CPU)**: Porting the code to FORTRAN and implementing with CPU (still serial solvers)

4. **Fortran + PETSc (MPI + CUDA)**: Using distributed memory parallelism and linear solvers enabled by GPU

Again, each stage solves the same mathematical problem to ensure consistency and traceability.



## Repository Structure

```text
poisson-fem/
│
├── README.md
├── python/
│   ├── 1d_cpu.ipynb
│   ├── 2d_cpu.ipynb
│   ├── cpu_vs_gpu_cupy.ipynb
│
├── fortran/
│   ├── 1D_scratch.f90
│   ├── 1D_random_grid.f90
│   ├── 2d_no_mpi.f90
│   ├── 2d_cpu_parallelize.F90
│   ├── petsc_gpu.f90
│   └── Makefile
│
├── scripts/
│   ├── run_python_cpu.sh
│   ├── run_python_gpu.sh
│   ├── run_petsc_cpu.sh
│   ├── run_petsc_gpu.sh
│
├── results/
│   ├── figures/
│   ├── timings/
│   └── petsc_logs/
│
└── environment/
    ├── requirements.txt
    ├── cuda_info.txt
    └── petsc_config.txt
```

