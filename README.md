# Finite Element Solution of the Poisson Equation  

## Overview

This repository documents the development of a **finite element solver for the Poisson equation**, starting from simple CPU-based implementations and extending to parallel and later GPU-accelerated solvers using PETSc.

The same mathematical problem is solved across multiple programming languages and hardware architectures in order to verify the numerical correctness and study performance scaling in MPI and GPU solvers.


## Problem Description

We consider the Poisson equation

$$\nabla^2 u(\mathbf{x}) = f(\mathbf{x})$$

on one- and two-dimensional domains with Dirichlet boundary conditions. The problem is discretized using the finite element method (FEM) with low-order basis functions and explicit matrix assembly. Full mathematical formulation and discretization details are provided in the documentation.



## Implementation Workflow

The solver is developed in the following manner (Wherever needed, job scripts have also been attached with the name corresponding to the source file replace with ".sh" in the directory: poisson-fem/jobscripts):

1. **Python (NumPy, CPU)**

2. **Python (CuPy, GPU)**: GPU Implementation in Python (Without MPI)

3. **Fortran + LAPACK (CPU)**: Porting the code to FORTRAN and implementing with CPU (still serial solvers). The first three files in the poisson-fem/fortran: "1D_scratch.f90", "1D_random_grid.f90" and "2d_no_mpi.f90" explores this progressively.

4. **Fortran + PETSc (MPI + CUDA)**: Using distributed memory parallelism and linear solvers enabled by GPU. In the directory "poisson-fem/fortran", we implement this progressively. Beginning with mpi paralleization introduction using CPU through the code in "2d_cpu_parallelize.F90" as the file name suggests. We then move on to the GPU implementation of the same code. Leveraging PETSc would be very useful here since it helps us skip many steps to directly implement the poisson solver giving the flexibility to execute it with mutiple CPU or GPU nodes. "petsc_gpu.F90" does the work effeciently.

   The _final opus_ of the repo ends with leveraging all the PETSc flags and features at the best and solving the poisson solver with the flexibility of execution in CPU or GPU **and having multiple source terms**. The addition of the **multiple soruce terms** instead of one makes the task bit more difficult but using additional variables, redeifing the communicators and making changes to the corresponding job scripts help achieve the work.

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
│   ├── petsc_gpu.F90
│   └── multi_probs.F90
│
├── scripts/
│   ├── 1D_scratch.sh
│   ├── 2d_cpu_parallelize.sh
│   ├── petsc_gpu.sh
│   ├── multi_probs.sh
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
## Motivation and Goal

Many finite element formulations produce large sparse linear systems that become computationally expensive as the problem size increases. While dense matrix implementations could be used for small-scale problems, scientific simulations works with sparse storage formats and parallel algorithms to efficiently use new-age computing hardware.

This project explores the complete workflow of transitioning from a basic serial finite element implementation to a high-performance solver by leveraging CuPy for sparse matrix assembly, Krylov-based linear solvers, and MPI-enabled distributed parallelism. Building upon this foundation, GPU acceleration is introduced through PETSc's CUDA backend, enabling faster execution on heterogeneous HPC systems without requiring significant changes to the application code.

The techniques developed here form the computational basis for accelerating plasma physics applications, particularly gyrokinetic simulations used in magnetic confinement fusion research. Codes such as G2C3 solve Poisson-like equations as part of their self-consistent field calculations, making sparse linear algebra and multi-GPU parallelism important for reducing simulation time and helping with higher-resolution studies of tokamak plasmas.

Although this repository demonstrates these ideas using the classical Poisson equation, the methodologies presented:
- sparse matrix assembly 
- PETSc-based linear solvers,
- MPI parallelization and
- GPU acceleration

are directly applicable to gyrokinetic Poisson solvers used in plasma physics and tokamak simulation codes.

**The ultimate goal is to extend these methodologies to the gyrokinetic Poisson solver used in tokamak plasma simulations. In such simulations, the computational domain can be decomposed into multiple independent subdomains (e.g., 32 or more), allowing the Poisson equation to be solved concurrently across distributed CPU cores or GPUs using MPI. This parallel decomposition significantly reduces computation time and forms a key component of scalable gyrokinetic plasma simulation codes for magnetic confinement fusion.**


