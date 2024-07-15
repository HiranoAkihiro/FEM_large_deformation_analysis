module mod_solver
    use mod_util
    implicit none
    
contains

    subroutine solver(mesh, var) !//TODO CG法を直接実装する
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        real(8), allocatable :: b(:)

        ! allocate(b(ndof*mesh%nnode))
        ! mat%dence%B = var%B
        ! b = var%B
        ! mat%dence%X = 0.0d0

        ! ! call Gauss_Jordan(mat%dence%A, b, mat%dence%X, ndof*mesh%nnode)
        ! call solver_CG_dense(mat%dence%A, mat%dence%B, mat%dence%X, ndof*mesh%nnode) !//[x] OK! solver_CG_dense
        ! var%X = mat%dence%X
        !call soild_debug_header("solver")

        call monolis_set_method(mat, monolis_iter_CG)
        call monolis_set_precond(mat, monolis_prec_DIAG)
        call monolis_set_maxiter(mat, 100000)
        call monolis_set_tolerance(mat, 1.0d-8)
        !call monolis_param_set_is_scaling(mat, .false.)
        !call monolis_param_set_is_reordering(mat, .false.)
        !call monolis_param_set_is_debug(mat, .true.)
        call monolis_show_timelog(mat, .true.)
        call monolis_show_iterlog(mat, .true.)
        call monolis_show_summary(mat, .true.)

        call monolis_solve_R(mat, com, var%B, var%X)
    !    call soild_plot_solver(mat%PRM%curiter, mat%PRM%curresid)

    !    if(mat%PRM%curresid > mat%PRM%tol)then
    !      if(mat%COM%myrank == 0) write(*,"(a)") "*** ERROR: monolis solver is not converge"
    !      stop
    !    endif
    end subroutine solver

    subroutine solver_CG_dense(A,b,x,n)
        implicit none
        integer(4) :: i,j,n
        real(8), intent(in) :: A(:,:) !　係数行列
        real(8), intent(in) :: b(:)   !　右辺ベクトル
        real(8), intent(inout) :: x(:)   !　解ベクトル
        real(8), allocatable :: r(:)   !　残差ベクトル
        real(8), allocatable :: p(:)   !　
        real(8), allocatable :: q(:)   !　
        real(8) :: bnrm2,rho,beta,rho1,dnrm2,c1,alpha,resid !　
        real(8) :: eps
        integer(4) :: iter,itermax
        character :: esc*1 = char(27)
        logical :: iter_log = .true.

        itermax=100000
        eps=1.0E-8
    
        x = 0.0d0
        allocate(r(n), source = 0.0d0)
        allocate(p(n), source = 0.0d0)
        allocate(q(n), source = 0.0d0)
    
        do i=1,n
            r(i) = b(i)
            do j=1,n
                r(i) = r(i) - A(i,j)*x(j)
            enddo
        enddo
    
        bnrm2 = 0.0d0
        do i=1,n
            bnrm2 = bnrm2 + b(i) ** 2.0d0
        enddo
    
        rho1 = 0.0d0
        do iter=1,itermax
    
            rho = 0.0d0
            do i=1,n
                rho = rho + r(i)*r(i)
            enddo
    
            if(iter.eq.1) then
                do i=1,n
                    p(i) = r(i)
                enddo
            else
                beta = rho / rho1
                do i=1,n
                    p(i) = r(i) + beta * p(i)
                enddo
            endif
    
            do i=1,n
                q(i) = 0.0d0
                do j=1,n
                    q(i) = q(i) + A(i,j) * p(j)
                enddo
            enddo
    
    
            c1=0.0d0
            do i=1,n
                c1 = c1 + p(i) * q(i)
            enddo
            alpha = rho / c1
    
            do i=1,n
                x(i) = x(i) + alpha * p(i)
                r(i) = r(i) - alpha * q(i)
            enddo
    
            dnrm2 = 0.0d0
            do i=1,n
                dnrm2 = dnrm2 + r(i) ** 2
            enddo
    
            resid = dsqrt(dnrm2/bnrm2)
            if(iter_log)write(*,'(i0,a,es20.6)')iter,' ',resid
    
            if(resid.le.eps) then
                write(*,'(a, i0)')'The number of solver iteration is ',iter
                exit
            endif
            rho1 =rho
    
            if(iter == itermax) then
                write(*,'(a,a)')esc//"[31m"//'<ITERATIVE SOLVER ERROR>'//esc//"[0m"&
                            ,' reached maximum iterations'
                stop
            endif
        enddo
    end subroutine solver_CG_dense

    subroutine Gauss_Jordan(A,b,x,n)
        implicit none
        integer(4) :: i,j,k,n
        real(8), intent(in) :: A(:,:) !　係数行列
        real(8), intent(in) :: b(:)   !　右辺ベクトル
        real(8), allocatable :: a_temp(:,:)
        real(8), allocatable :: b_temp(:)
        real(8), intent(inout) :: x(:)   !　解ベクトル
        real(8) :: multiplier

        allocate(a_temp(n,n), source = 0.0d0)
        allocate(b_temp(n), source = 0.0d0)

        a_temp = A
        b_temp = b

        do i=1,n-1
            do j=i+1,n
                multiplier = a_temp(j,i)/a_temp(i,i)
                do k=i,n
                    a_temp(j,k) = a_temp(j,k)-multiplier*a_temp(i,k)
                enddo
                b_temp(j) = b_temp(j)-multiplier*b(i)
            enddo
        enddo

        do i=n,1,-1
            x(i) = b_temp(i)/a_temp(i,i)
            do j=i+1,n
                x(i) = x(i) - a_temp(i,j)*x(j)/a_temp(i,i)
            end do
        end do
    end subroutine Gauss_Jordan

    function is_convergence(mesh, var, step)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        integer(4) :: step, i
        real(8), save :: b0nrm
        real(8) :: bnrm, rnrm, rnrmmax, qnrm
        logical :: is_convergence
        character :: esc*1 = char(27)

        is_convergence = .false.

        bnrm = 0.0d0
        rnrm = 0.0d0
        qnrm = 0.0d0
        bnrm = 0.0d0

        do i = 1, ndof*mesh%nnode
            qnrm = qnrm + var%q(i)*var%q(i)
            bnrm = bnrm + var%b(i)*var%b(i)
        enddo

        if(step == 1)then
            b0nrm = bnrm
        write(*,"(a,1pe12.5)")esc//"[36m"//"  ** NR        b0nrm: "//esc//"[0m", dsqrt(b0nrm)
        else
            !rnrm    = dsqrt(bnrm/b0nrm)
            rnrm    = dsqrt(bnrm/qnrm)
            rnrmmax = dabs(maxval(var%B))
            write(*,"(a,1pe12.5,a,1pe12.5)")esc//"[36m"//"  ** NR     residual: "//esc//"[0m", rnrm, ", ", rnrmmax
            if(rnrm < 1.0d-6 .or. rnrmmax < 1.0d-8) is_convergence = .true.
        endif
    end function is_convergence
end module mod_solver