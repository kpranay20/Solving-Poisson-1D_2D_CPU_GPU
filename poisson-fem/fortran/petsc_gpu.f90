#include "petsc/finclude/petscksp.h"
program petsc_gpu
  use petscksp
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: num_grids = 150
  integer, parameter :: num_points = (num_grids + 1)**2
  integer, parameter :: num_triangles = 2 * num_grids**2

  real(dp) :: points(num_points,2)
  integer :: triangles(num_triangles,3)
  logical :: boundary_nodes(num_points)
  integer :: i,j,t,tri_count
  integer :: global_i, global_j
  real(dp) :: x,y,dx,dy
  real(dp) :: coords(3,2), grads(3,2), area
  real(dp) :: kt(3,3), ft(3), f_values(3)
  integer :: point_indices(num_grids+1, num_grids+1)
  integer :: bl, br, tr, tl
  real(dp) :: start_time, end_time
  real(dp) :: norm_err, norm_u, norm_u_exact, err
  real(dp), allocatable :: u(:), u_exact(:)

  Mat :: A
  Vec :: b, x_vec
  KSP :: ksp
  PetscErrorCode :: ierr
  PetscInt :: n

  call cpu_time(start_time)
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  n = num_points
  allocate(u(n), u_exact(n))

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
        triangles(tri_count,:) = [bl, br, tr]
        tri_count = tri_count + 1
        triangles(tri_count,:) = [bl, tr, tl]
        tri_count = tri_count + 1
     end do
  end do

  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  call MatSetType(A, MATSEQAIJCUSPARSE, ierr) ! GPU matrix
  call MatSetUp(A, ierr)

  call VecCreate(PETSC_COMM_WORLD, b, ierr)
  call VecSetSizes(b, PETSC_DECIDE, n, ierr)
  call VecSetType(b, VECCUDA, ierr)           
  call VecDuplicate(b, x_vec, ierr)          

  do t = 1, num_triangles
     do i = 1, 3
        coords(i,:) = points(triangles(t,i),:)
     end do
     call compute_basis_functions(coords, grads, area)
     kt = matmul(grads, transpose(grads)) * area

     do i = 1,3
        x = coords(i,1)
        y = coords(i,2)
        f_values(i) = source_term(x,y)
     end do
     ft = f_values * (area/3.0_dp)

     do i = 1,3
        global_i = triangles(t,i)
        call VecSetValue(b, global_i-1, ft(i), ADD_VALUES, ierr)
        do j = 1,3
           global_j = triangles(t,j)
           call MatSetValue(A, global_i-1, global_j-1, kt(i,j), ADD_VALUES, ierr)
        end do
     end do
  end do

  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  call VecAssemblyBegin(b, ierr)
  call VecAssemblyEnd(b, ierr)

  boundary_nodes = .false.
  do i = 1, n
     x = points(i,1)
     y = points(i,2)
     if (abs(x)<1.0e-12_dp .or. abs(x-1.0_dp)<1.0e-12_dp .or. &
         abs(y)<1.0e-12_dp .or. abs(y-1.0_dp)<1.0e-12_dp) then
        boundary_nodes(i) = .true.
     end if
  end do

  do i = 1, n
     if (boundary_nodes(i)) then
        call MatZeroRows(A,1,[i-1],1.0_dp,x_vec,b,ierr)
     end if
  end do

  call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call KSPSolve(ksp, b, x_vec, ierr)

  call VecGetValues(x_vec, n, [(i-1,i=1,n)], u, ierr)

  do i = 1, n
     x = points(i,1)
     y = points(i,2)
     u_exact(i) = analytical_solution(x,y)
  end do

  norm_u_exact = sqrt(sum(u_exact**2))
  norm_err = sqrt(sum((u - u_exact)**2))/norm_u_exact
  print '(A,F10.4)', 'L2 Relative Error: ', norm_err

  call KSPDestroy(ksp, ierr)
  call VecDestroy(x_vec, ierr)
  call VecDestroy(b, ierr)
  call MatDestroy(A, ierr)
  call PetscFinalize(ierr)

  call cpu_time(end_time)
  print '(A,F8.2,A)','GPU computation time: ',end_time-start_time,' seconds'

contains

  function source_term(x,y) result(f)
    real(dp), intent(in)::x,y
    real(dp)::f
    real(dp), parameter::pi=3.141592653589793_dp
    f = sin(pi*x)*sin(pi*y)
  end function source_term

  function analytical_solution(x,y) result(u_anal)
    real(dp), intent(in)::x,y
    real(dp)::u_anal
    real(dp), parameter::pi=3.141592653589793_dp
    u_anal = (1.0_dp/(2.0_dp*pi**2))*sin(pi*x)*sin(pi*y)
  end function analytical_solution

  subroutine compute_basis_functions(coords,grads,area)
    real(dp), intent(in)::coords(3,2)
    real(dp), intent(out)::grads(3,2),area
    real(dp)::x1,x2,x3,y1,y2,y3,detJ
    x1=coords(1,1); y1=coords(1,2)
    x2=coords(2,1); y2=coords(2,2)
    x3=coords(3,1); y3=coords(3,2)
    detJ=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
    area=0.5_dp*abs(detJ)
    grads(1,:)=[(y2-y3),(x3-x2)]/detJ
    grads(2,:)=[(y3-y1),(x1-x3)]/detJ
    grads(3,:)=[(y1-y2),(x2-x1)]/detJ
  end subroutine compute_basis_functions

end program petsc_gpu
