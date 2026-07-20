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
## Motivation

Many finite element formulations produce large sparse linear systems that become computationally expensive as the problem size increases. While dense matrix implementations could be used for small-scale problems, scientific simulations works with sparse storage formats and parallel algorithms to efficiently use new-age computing hardware.

This project explores the complete workflow of transitioning from a basic serial finite element implementation to a high-performance solver by leveraging CuPy for sparse matrix assembly, Krylov-based linear solvers, and MPI-enabled distributed parallelism. Building upon this foundation, GPU acceleration is introduced through PETSc's CUDA backend, enabling faster execution on heterogeneous HPC systems without requiring significant changes to the application code.

The techniques developed here form the computational basis for accelerating plasma physics applications, particularly gyrokinetic simulations used in magnetic confinement fusion research. Codes such as G2C3 solve Poisson-like equations as part of their self-consistent field calculations, making sparse linear algebra and multi-GPU parallelism important for reducing simulation time and helping with higher-resolution studies of tokamak plasmas.

Although this repository demonstrates these ideas using the classical Poisson equation, the methodologies presented:
- sparse matrix assembly 
- PETSc-based linear solvers,
- MPI parallelization and
- GPU acceleration

are directly applicable to gyrokinetic Poisson solvers used in plasma physics and tokamak simulation codes.
