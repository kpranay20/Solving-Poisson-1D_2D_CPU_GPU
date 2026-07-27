#include "petsc/finclude/petscksp.h"
program petsc_parallel
  use petsc
  use petscksp
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: num_grids = 150
  integer, parameter :: num_points = (num_grids + 1)**2
  integer, parameter :: num_triangles = 2 * num_grids**2

  real(dp) :: points(num_points, 2)
  integer :: triangles(num_triangles, 3)
  logical :: boundary_nodes(num_points)
  integer :: i, j, t, tri_count
  integer :: global_i, global_j
  real(dp) :: x, y, dx, dy
  real(dp) :: coords(3, 2), grads(3, 2), area
  real(dp) :: kt(3, 3), ft(3), f_values(3)
  integer :: point_indices(num_grids + 1, num_grids + 1)
  integer :: bl, br, tr, tl
  real(dp) :: start_time, end_time
  real(dp) :: norm_err
  PetscReal :: norm_u_exact_p, norm_diff_p
  integer :: rank, nprocs
  integer :: nb, k

  Mat :: A
  Vec :: b, x_vec, u_exact_vec, diff_vec
  KSP :: ksp
  IS  :: is_bnd
  PetscErrorCode :: ierr
  PetscInt :: n, istart, iend
  PetscInt, allocatable :: bnd_idx(:)

  Vec :: x_seq
  VecScatter :: scat

  integer :: fh, gi
  PetscScalar, pointer :: xarr(:)
  real(dp) :: ue, err
  integer :: gi_p, gj_p    

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, nprocs, ierr)

  n = num_points

  dx = 1.0_dp / real(num_grids, dp)
  dy = 1.0_dp / real(num_grids, dp)

  tri_count = 1
  do j = 0, num_grids
    do i = 0, num_grids
      points(j*(num_grids+1)+i+1,1) = real(i,dp)*dx
      points(j*(num_grids+1)+i+1,2) = real(j,dp)*dy
      point_indices(j+1,i+1) = j*(num_grids+1)+i+1
    end do
  end do

  do i = 1, num_grids
    do j = 1, num_grids
      bl = point_indices(j,i)
      br = point_indices(j,i+1)
      tr = point_indices(j+1,i+1)
      tl = point_indices(j+1,i)

      triangles(tri_count,:) = [bl, br, tr]; tri_count = tri_count + 1
      triangles(tri_count,:) = [bl, tr, tl]; tri_count = tri_count + 1
    end do
  end do

  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  call MatSetType(A, MATAIJ, ierr)   
  call MatSetUp(A, ierr)

  call VecCreate(PETSC_COMM_WORLD, b, ierr)
  call VecSetSizes(b, PETSC_DECIDE, n, ierr)
  call VecSetFromOptions(b, ierr)

  call VecDuplicate(b, x_vec, ierr)
  call VecDuplicate(b, u_exact_vec, ierr)
  call VecDuplicate(b, diff_vec, ierr)

  do t = 1, num_triangles
    do i = 1, 3
      coords(i, :) = points(triangles(t, i), :)
    end do

    call compute_basis_functions(coords, grads, area)
    kt = matmul(grads, transpose(grads)) * area

    do i = 1, 3
      x = coords(i,1)
      y = coords(i,2)
      f_values(i) = source_term(x, y)
    end do
    ft = f_values * (area / 3.0_dp)

    do i = 1, 3
      global_i = triangles(t, i)
      gi_p = global_i - 1
      call VecSetValue(b, gi_p, ft(i), ADD_VALUES, ierr)
      do j = 1, 3
        global_j = triangles(t, j)
        gj_p = global_j - 1
        call MatSetValue(A, gi_p, gj_p, kt(i,j), ADD_VALUES, ierr)
      end do
    end do
  end do

  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  call VecAssemblyBegin(b, ierr)
  call VecAssemblyEnd(b, ierr)

  boundary_nodes = .false.
  do i = 1, num_points
    x = points(i,1); y = points(i,2)
    if (abs(x) < 1.0e-12_dp .or. abs(x-1.0_dp) < 1.0e-12_dp .or. &
        abs(y) < 1.0e-12_dp .or. abs(y-1.0_dp) < 1.0e-12_dp) then
      boundary_nodes(i) = .true.
    end if
  end do

  nb = count(boundary_nodes)
  allocate(bnd_idx(nb))
  k = 0
  do i = 1, num_points
    if (boundary_nodes(i)) then
      k = k + 1
      bnd_idx(k) = i-1   ! PETSc being 0-based
    end if
  end do

  call ISCreateGeneral(PETSC_COMM_WORLD, nb, bnd_idx, PETSC_COPY_VALUES, is_bnd, ierr)
  call MatZeroRowsIS(A, is_bnd, 1.0d0, x_vec, b, ierr)
  call ISDestroy(is_bnd, ierr)
  deallocate(bnd_idx)

  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call cpu_time(start_time)

  call KSPSolve(ksp, b, x_vec, ierr)
  call cpu_time(end_time)

  call VecGetOwnershipRange(u_exact_vec, istart, iend, ierr)  ! [istart, iend) 0-based
  do i = istart, iend-1
    x = points(i+1,1)  ! convert back to 1-based for our arrays
    y = points(i+1,2)
    call VecSetValue(u_exact_vec, i, analytical_solution(x,y), INSERT_VALUES, ierr)
  end do
  call VecAssemblyBegin(u_exact_vec, ierr)
  call VecAssemblyEnd(u_exact_vec, ierr)

  call VecWAXPY(diff_vec, -1.0d0, u_exact_vec, x_vec, ierr)  ! diff = x + (-1) * u_exact
  call VecNorm(diff_vec, NORM_2, norm_diff_p, ierr)
  call VecNorm(u_exact_vec, NORM_2, norm_u_exact_p, ierr)
  if (norm_u_exact_p > 0.0d0) then
    norm_err = real(norm_diff_p / norm_u_exact_p, dp)
  else
    norm_err = 0.0_dp
  end if

  if (rank == 0) then
    write(*,'(A,E12.4)') 'L2 Relative Error: ', norm_err
  end if

  call VecScatterCreateToZero(x_vec, scat, x_seq, ierr)
  call VecScatterBegin(scat, x_vec, x_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd  (scat, x_vec, x_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)

  if (rank == 0) then
    call VecGetArrayF90(x_seq, xarr, ierr)

    open(unit=10, file='2D_cpu_parallelize.csv', status='replace', action='write')
    write(10,'(A)') 'Point_Index,X,Y,U_Numerical,U_Analytical,Error'
    do gi = 1, n
      x = points(gi,1); y = points(gi,2)
      ue = analytical_solution(x,y)
      err = abs(real(xarr(gi),dp) - ue)
      write(10,'(I0,5(",",ES15.8))') gi, x, y, real(xarr(gi),dp), ue, err
    end do
    close(10)

    call VecRestoreArrayF90(x_seq, xarr, ierr)
  end if

  call VecScatterDestroy(scat, ierr)
  call VecDestroy(x_seq, ierr)

  ! Clean
  call KSPDestroy(ksp, ierr)
  call VecDestroy(diff_vec, ierr)
  call VecDestroy(u_exact_vec, ierr)
  call VecDestroy(x_vec, ierr)
  call VecDestroy(b, ierr)
  call MatDestroy(A, ierr)

  call PetscFinalize(ierr)
  if (rank == 0) then
    print '(A,F8.2,A)', 'CPU computation time: ', end_time - start_time, ' seconds'
  end if

contains

  function source_term(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    real(dp), parameter :: pi = 3.141592653589793_dp
    f = sin(pi*x)*sin(pi*y)
  end function source_term

  function analytical_solution(x, y) result(u_anal)
    real(dp), intent(in) :: x, y
    real(dp) :: u_anal
    real(dp), parameter :: pi = 3.141592653589793_dp
    u_anal = (1.0_dp/(2.0_dp*pi**2))*sin(pi*x)*sin(pi*y)
  end function analytical_solution

  subroutine compute_basis_functions(coords, grads, area)
    real(dp), intent(in) :: coords(3,2)
    real(dp), intent(out) :: grads(3,2), area
    real(dp) :: x1,x2,x3,y1,y2,y3
    real(dp) :: detJ

    x1 = coords(1,1); y1 = coords(1,2)
    x2 = coords(2,1); y2 = coords(2,2)
    x3 = coords(3,1); y3 = coords(3,2)

    detJ = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
    area = 0.5_dp * abs(detJ)

    grads(1,:) = [(y2 - y3), (x3 - x2)] / detJ
    grads(2,:) = [(y3 - y1), (x1 - x3)] / detJ
    grads(3,:) = [(y1 - y2), (x2 - x1)] / detJ
  end subroutine compute_basis_functions

end program petsc_parallel
