#include "petsc/finclude/petscksp.h"
program multi_gpu
  use petsc
  use petscksp
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: num_problems = 2   ! change ONLY this when adding problems
  integer, parameter :: num_grids = 1000

  integer, parameter :: num_points = (num_grids + 1)**2
  integer, parameter :: num_triangles = 2 * num_grids**2

  !MPI
  integer :: rank, nprocs
  integer :: color, subRank, subSize
  MPI_Comm :: subComm

  !MESH
  real(dp) :: points(num_points,2)
  integer :: triangles(num_triangles,3)
  integer :: point_indices(num_grids+1, num_grids+1)
  logical :: boundary_nodes(num_points)

  real(dp) :: coords(3,2), grads(3,2), area
  real(dp) :: kt(3,3), ft(3), f_values(3)

  !VARIABLES
  integer :: i,j,t,tri_count
  integer :: global_i, global_j
  real(dp) :: x,y,dx,dy
  real(dp) :: ue, err

  !PETSc
  Mat :: A
  Vec :: b, x_vec, u_exact_vec, diff_vec
  KSP :: ksp
  IS  :: is_bnd

  PetscViewer :: viewer
  character(len=50) :: logname
  PetscErrorCode :: ierr
  PetscInt :: n

  !INIT
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  call PetscLogDefaultBegin(ierr)

  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, nprocs, ierr)

  !SPLIT
  color = rank / (nprocs / num_problems)
  call MPI_Comm_split(PETSC_COMM_WORLD, color, rank, subComm, ierr)

  call MPI_Comm_rank(subComm, subRank, ierr)
  call MPI_Comm_size(subComm, subSize, ierr)

  if (subRank == 0) then
    print *, "Problem", color, "running with", subSize, "ranks"
  end if

  !GRID
  n = num_points
  dx = 1.0_dp / real(num_grids, dp)
  dy = dx

  tri_count = 1
  do j = 0, num_grids
    do i = 0, num_grids
      points(j*(num_grids+1)+i+1,:) = [real(i,dp)*dx, real(j,dp)*dy]
      point_indices(j+1,i+1) = j*(num_grids+1)+i+1
    end do
  end do

  do i = 1, num_grids
    do j = 1, num_grids
      triangles(tri_count,:) = [point_indices(j,i), point_indices(j,i+1), point_indices(j+1,i+1)]
      tri_count = tri_count + 1
      triangles(tri_count,:) = [point_indices(j,i), point_indices(j+1,i+1), point_indices(j+1,i)]
      tri_count = tri_count + 1
    end do
  end do

  !PETSc SETUP
  call MatCreate(subComm, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  call MatSetFromOptions(A, ierr)
  call MatSetUp(A, ierr)

  call VecCreate(subComm, b, ierr)
  call VecSetSizes(b, PETSC_DECIDE, n, ierr)
  call VecSetFromOptions(b, ierr)

  call VecDuplicate(b, x_vec, ierr)
  call VecDuplicate(b, u_exact_vec, ierr)
  call VecDuplicate(b, diff_vec, ierr)
  !ASSEMBLY
  do t = 1, num_triangles

    do i = 1,3
      coords(i,:) = points(triangles(t,i),:)
    end do

    call compute_basis(coords, grads, area)
    kt = matmul(grads, transpose(grads)) * area

    do i = 1,3
      x = coords(i,1)
      y = coords(i,2)
      f_values(i) = source_term(x,y,color)
    end do

    ft = f_values * (area/3.0_dp)

    do i = 1,3
      global_i = triangles(t,i) - 1
      call VecSetValue(b, global_i, ft(i), ADD_VALUES, ierr)

      do j = 1,3
        global_j = triangles(t,j) - 1
        call MatSetValue(A, global_i, global_j, kt(i,j), ADD_VALUES, ierr)
      end do
    end do

  end do

  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  call VecAssemblyBegin(b, ierr)
  call VecAssemblyEnd(b, ierr)

  !BC
  boundary_nodes = .false.
  do i = 1, num_points
    x = points(i,1); y = points(i,2)
    if (x==0 .or. x==1 .or. y==0 .or. y==1) boundary_nodes(i)=.true.
  end do

  call apply_bc(A, b, boundary_nodes, subComm)

  !SOLVE
  call KSPCreate(subComm, ksp, ierr)
  call KSPSetOperators(ksp, A, A, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call KSPSolve(ksp, b, x_vec, ierr)

 !Log Files
  write(logname,'(A,I0,A)') 'log_problem_', color, '.txt'

  call PetscViewerASCIIOpen(subComm, logname, viewer, ierr)
  call PetscLogView(viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  !ERROR
  call compute_error(points, x_vec, u_exact_vec, diff_vec, color, subComm)

  !OUTPUT
  call write_output(points, x_vec, color, subComm)

  !CLEANUP
  call KSPDestroy(ksp, ierr)
  call VecDestroy(diff_vec, ierr)
  call VecDestroy(u_exact_vec, ierr)
  call VecDestroy(x_vec, ierr)
  call VecDestroy(b, ierr)
  call MatDestroy(A, ierr)

  call MPI_Comm_free(subComm, ierr)
  call PetscFinalize(ierr)

contains

!SOURCE
function source_term(x,y,prob) result(f)
  real(dp), intent(in) :: x,y
  integer, intent(in) :: prob
  real(dp) :: f
  real(dp), parameter :: pi = 3.141592653589793_dp

  select case(prob)
  case(0)
    f = sin(pi*x)*sin(pi*y)
  case(1)
    f = 2*pi**2*(sin(pi*x)*sin(pi*y)+sin(2*pi*x)*sin(2*pi*y))
  case default
    f = 0.0_dp
  end select
end function

!ANALYTICAL
function analytical_solution(x,y,prob) result(u)
  real(dp), intent(in) :: x,y
  integer, intent(in) :: prob
  real(dp) :: u
  real(dp), parameter :: pi = 3.141592653589793_dp

  select case(prob)
  case(0)
    u = sin(pi*x)*sin(pi*y)/(2*pi**2)
  case(1)
    u = sin(pi*x)*sin(pi*y)+0.25*sin(2*pi*x)*sin(2*pi*y)
  case default
    u = 0.0_dp
  end select
end function

!BASIS
subroutine compute_basis(coords, grads, area)
  real(dp), intent(in) :: coords(3,2)
  real(dp), intent(out) :: grads(3,2), area
  real(dp) :: x1,x2,x3,y1,y2,y3, detJ

  x1=coords(1,1); y1=coords(1,2)
  x2=coords(2,1); y2=coords(2,2)
  x3=coords(3,1); y3=coords(3,2)

  detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
  area = 0.5_dp*abs(detJ)

  grads(1,:) = [(y2-y3),(x3-x2)]/detJ
  grads(2,:) = [(y3-y1),(x1-x3)]/detJ
  grads(3,:) = [(y1-y2),(x2-x1)]/detJ
end subroutine

!BC
subroutine apply_bc(A,b,boundary,comm)
  Mat :: A
  Vec :: b
  logical :: boundary(:)
  MPI_Comm :: comm
  integer :: i, nb
  PetscInt, allocatable :: idx(:)
  IS :: is_bnd
  PetscErrorCode :: ierr

  nb = count(boundary)
  allocate(idx(nb))
  nb = 0
  do i=1,size(boundary)
    if (boundary(i)) then
      nb=nb+1
      idx(nb)=i-1
    end if
  end do

  call ISCreateGeneral(comm, nb, idx, PETSC_COPY_VALUES, is_bnd, ierr)
  call MatZeroRowsIS(A, is_bnd, 1.0d0, PETSC_NULL_VEC, b, ierr)
  call ISDestroy(is_bnd, ierr)
  deallocate(idx)
end subroutine

!ERROR
subroutine compute_error(points,x,u_exact,diff,prob,comm)
  real(dp) :: points(:,:)
  Vec :: x,u_exact,diff
  integer :: prob
  MPI_Comm :: comm
  PetscErrorCode :: ierr
  PetscInt :: i, istart,iend
  PetscReal :: nrm1,nrm2
  real(dp) :: xx,yy

  call VecGetOwnershipRange(x,istart,iend,ierr)
  do i=istart,iend-1
    xx=points(i+1,1); yy=points(i+1,2)
    call VecSetValue(u_exact,i,analytical_solution(xx,yy,prob),INSERT_VALUES,ierr)
  end do

  call VecAssemblyBegin(u_exact,ierr)
  call VecAssemblyEnd(u_exact,ierr)

  call VecWAXPY(diff,-1.0d0,u_exact,x,ierr)
  call VecNorm(diff,NORM_2,nrm1,ierr)
  call VecNorm(u_exact,NORM_2,nrm2,ierr)

  if (nrm2>0) then
    if (istart==0) print *, "Relative Error:", nrm1/nrm2
  end if
end subroutine

!OUTPUT
subroutine write_output(points,x,prob,comm)
  real(dp) :: points(:,:)
  Vec :: x
  integer :: prob
  MPI_Comm :: comm
  PetscErrorCode :: ierr
  Vec :: seq
  VecScatter :: scat
  PetscScalar, pointer :: arr(:)
  character(len=50) :: fname

  call VecScatterCreateToZero(x,scat,seq,ierr)
  call VecScatterBegin(scat,x,seq,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(scat,x,seq,INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(seq,arr,ierr)

  write(fname,'(A,I0,A)') 'results_problem_',prob,'.csv'
  open(10,file=fname,status='replace')

  do i=1,size(points,1)
    write(10,'(I0,3(",",ES15.8))') i,points(i,1),points(i,2),arr(i)
  end do

  close(10)
  call VecRestoreArrayF90(seq,arr,ierr)
end subroutine

end program
