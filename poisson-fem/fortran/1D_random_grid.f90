program fem_poisson
    implicit none
    integer, parameter :: n = 101   ! 'parameter' make 'n' un-assignable for later   
    integer :: i
    real(8) :: h(n), random_value
    real(8) :: x(n), ke(2,2), fe(2)
    real(8), dimension(n,n) :: K
    real(8), dimension(n)   :: F, u, u_exact


    ! node points
    do i = 1, n
        x(i) = (i-1) * h
    end do

    ! Local stiffness and load matrices
    ke = reshape([1d0, -1d0, -1d0, 1d0], shape(ke)) / h
    fe = - (h/2d0) * (/ 1d0, 1d0 /)

    ! Initialize
    K = 0.0d0
    F = 0.0d0

    ! Assembly of global matrices
    do i = 1, n-1
        K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + ke
        F(i:i+1) = F(i:i+1) + fe
    end do

    ! bound. condn: u(0)=0, u(1)=0
    K(1,:) = 0.0d0
    K(1,1) = 1.0d0
    F(1)   = 0.0d0

    K(n,:) = 0.0d0
    K(n,n) = 1.0d0
    F(n)   = 0.0d0

    ! Solve linear system Ku = F (Gauss elimination)
    call gauss_solve(K, F, u, n)

    ! Analytical solution
    do i = 1, n
        u_exact(i) = -0.5d0 * x(i) * (1.0d0 - x(i))
    end do

    ! Write results to file
    open(unit=10, file="solution.dat", status="replace", action="write")
    write(10, '(A)') "#   x        FEM_Solution     Analytical_Solution"
    do i = 1, n
        write(10,'(F10.5, 3X, F15.8, 3X, F15.8)') x(i), u(i), u_exact(i)
    end do
    close(10)

    print *, "Results written to solution.dat"

contains

    subroutine gauss_solve(A, b, x, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(inout) :: A(n,n)
        real(8), intent(inout) :: b(n)
        real(8), intent(out) :: x(n)
        integer :: i, j, k, maxrow
        real(8) :: tmp, factor

        ! Forward elimination
        do k = 1, n-1
            ! Pivoting: find row with largest absolute value in column k
            maxrow = k
            do i = k+1, n
                if (abs(A(i,k)) > abs(A(maxrow,k))) maxrow = i
            end do

            ! Swap rows if needed
            if (maxrow /= k) then
                A([k,maxrow],:) = A([maxrow,k],:)   ! swap rows in A
                tmp = b(k); b(k) = b(maxrow); b(maxrow) = tmp
            end if

            ! Elimination
            do i = k+1, n
                if (A(k,k) /= 0.0d0) then
                    factor = A(i,k) / A(k,k)
                    A(i,k:n) = A(i,k:n) - factor * A(k,k:n)
                    b(i) = b(i) - factor * b(k)
                end if
            end do
        end do

        ! Back substitution
        do i = n, 1, -1
            tmp = b(i)
            do j = i+1, n
                tmp = tmp - A(i,j) * x(j)
            end do
            if (A(i,i) /= 0.0d0) then
                x(i) = tmp / A(i,i)
            else
                x(i) = 0.0d0
            end if
        end do
    end subroutine gauss_solve

end program fem_poisson
