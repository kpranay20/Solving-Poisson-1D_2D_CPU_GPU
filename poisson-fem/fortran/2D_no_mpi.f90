!Run with:  gfortran 2D_no_mpi.f90 -llapack -lblas -o fem_poisson

program fem_poisson
    implicit none
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: num_grids = 50
    integer, parameter :: num_points = (num_grids + 1)**2
    integer, parameter :: num_triangles = 2 * num_grids**2
    
    real(dp) :: points(num_points, 2)
    integer :: triangles(num_triangles, 3)
    real(dp) :: K(num_points, num_points)
    real(dp) :: F(num_points)
    real(dp) :: u(num_points)
    real(dp) :: u_exact(num_points)
    logical :: boundary_nodes(num_points)
    
    integer :: i, j, idx, t, node
    integer :: global_i, global_j
    real(dp) :: x, y, dx, dy
    real(dp) :: coords(3, 2)
    real(dp) :: A(3, 3), rhs(3), sol(3)
    real(dp) :: grads(3, 2), area
    real(dp) :: kt(3, 3), ft(3)
    real(dp) :: f_values(3)
    real(dp) :: start_time, end_time
    real(dp) :: error_l2, norm_u, norm_u_exact
    integer :: point_indices(num_grids + 1, num_grids + 1)
    integer :: bl, br, tr, tl
    integer :: tri_count
    
    call cpu_time(start_time)
    
    dx = 1.0_dp / real(num_grids, dp)
    dy = 1.0_dp / real(num_grids, dp)
    
    idx = 1
    do j = 0, num_grids
        do i = 0, num_grids
            points(idx, 1) = real(i, dp) * dx
            points(idx, 2) = real(j, dp) * dy
            point_indices(j + 1, i + 1) = idx
            idx = idx + 1
        end do
    end do
    
    tri_count = 1
    do i = 1, num_grids
        do j = 1, num_grids
            bl = point_indices(j, i)
            br = point_indices(j, i + 1)
            tr = point_indices(j + 1, i + 1)
            tl = point_indices(j + 1, i)
            
            triangles(tri_count, 1) = bl
            triangles(tri_count, 2) = br
            triangles(tri_count, 3) = tr
            tri_count = tri_count + 1
            
            triangles(tri_count, 1) = bl
            triangles(tri_count, 2) = tr
            triangles(tri_count, 3) = tl
            tri_count = tri_count + 1
        end do
    end do
    
    K = 0.0_dp
    F = 0.0_dp
    
    do t = 1, num_triangles
        do i = 1, 3
            coords(i, 1) = points(triangles(t, i), 1)
            coords(i, 2) = points(triangles(t, i), 2)
        end do
        
        call compute_basis_functions(coords, grads, area)
        
        kt = matmul(grads, transpose(grads)) * area
        
        do i = 1, 3
            x = coords(i, 1)
            y = coords(i, 2)
            f_values(i) = source_term(x, y)
        end do
        ft = f_values * (area / 3.0_dp)
        
        do i = 1, 3
            global_i = triangles(t, i)
            F(global_i) = F(global_i) + ft(i)
            do j = 1, 3
                global_j = triangles(t, j)
                K(global_i, global_j) = K(global_i, global_j) + kt(i, j)
            end do
        end do
    end do
    
    boundary_nodes = .false.
    do i = 1, num_points
        x = points(i, 1)
        y = points(i, 2)
        if (abs(x) < 1.0e-12_dp .or. abs(x - 1.0_dp) < 1.0e-12_dp .or. &
            abs(y) < 1.0e-12_dp .or. abs(y - 1.0_dp) < 1.0e-12_dp) then
            boundary_nodes(i) = .true.
        end if
    end do
    
    do i = 1, num_points
        if (boundary_nodes(i)) then
            K(i, :) = 0.0_dp
            K(:, i) = 0.0_dp
            K(i, i) = 1.0_dp
            F(i) = 0.0_dp
        end if
    end do
    
    call solve_linear_system(K, F, u, num_points)
    
    call cpu_time(end_time)
    print '(A,F8.2,A)', 'CPU computation time: ', end_time - start_time, ' seconds'
    
    do i = 1, num_points
        x = points(i, 1)
        y = points(i, 2)
        u_exact(i) = analytical_solution(x, y)
    end do
    
    norm_u = sqrt(sum(u**2))
    norm_u_exact = sqrt(sum(u_exact**2))
    error_l2 = sqrt(sum((u - u_exact)**2)) / norm_u_exact
    print '(A,F8.4)', 'L2 Relative Error: ', error_l2
    
    call write_results_csv(points, u, u_exact, num_points)
    
    print *, 'Results written to 2D_no_mpi.csv'

contains

    function source_term(x, y) result(f)
        real(dp), intent(in) :: x, y
        real(dp) :: f
        real(dp), parameter :: pi = 3.141592653589793_dp
        f = sin(pi * x) * sin(pi * y)
    end function source_term
    
    function analytical_solution(x, y) result(u_anal)
        real(dp), intent(in) :: x, y
        real(dp) :: u_anal
        real(dp), parameter :: pi = 3.141592653589793_dp
        u_anal = (1.0_dp / (2.0_dp * pi**2)) * sin(pi * x) * sin(pi * y)
    end function analytical_solution
    
    subroutine compute_basis_functions(coords, grads, area)
        real(dp), intent(in) :: coords(3, 2)
        real(dp), intent(out) :: grads(3, 2), area
        real(dp) :: A(3, 3), rhs(3), sol(3)
        integer :: k, info
        integer :: ipiv(3)
        
        A(:, 1) = 1.0_dp
        A(:, 2:3) = coords
        
        do k = 1, 3
            rhs = 0.0_dp
            rhs(k) = 1.0_dp
            sol = rhs
            
            call dgesv(3, 1, A, 3, ipiv, sol, 3, info)
            if (info /= 0) then
                print *, 'Error in linear solve for basis functions'
                stop
            end if
            
            grads(k, :) = sol(2:3)
            
            ! Reset A for next iteration
            A(:, 1) = 1.0_dp
            A(:, 2:3) = coords
        end do
        
        area = 0.5_dp * abs(A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
                            A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
                            A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)))
    end subroutine compute_basis_functions
    
    !LAPACK
    subroutine solve_linear_system(A, b, x, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n, n)
        real(dp), intent(inout) :: b(n)
        real(dp), intent(out) :: x(n)
        integer :: ipiv(n), info
        
        x = b
        call dgesv(n, 1, A, n, ipiv, x, n, info)
        
        if (info /= 0) then
            print *, 'Error in linear system solve, info = ', info
            stop
        end if
    end subroutine solve_linear_system
    
    subroutine write_results_csv(points, u, u_exact, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: points(n, 2), u(n), u_exact(n)
        integer :: i
        
        open(unit=10, file='2D_no_mpi.csv', status='replace')
        write(10, '(A)') 'Point_Index,X_Coordinate,Y_Coordinate,U_Numerical,U_Analytical,Error'
        
        do i = 1, n
            write(10, '(I0,5(",",ES15.8))') i, points(i, 1), points(i, 2), &
                                           u(i), u_exact(i), abs(u(i) - u_exact(i))
        end do
        
        close(10)
    end subroutine write_results_csv

end program fem_poisson