module mod_element
    use mod_util
    use mod_debug
    implicit none

    ! private

    real(8), parameter :: gsp(3,8) = reshape([ &
    -0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0,-0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0,-0.577350269189626d0, &
    -0.577350269189626d0,-0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0,-0.577350269189626d0, 0.577350269189626d0, &
    -0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0, &
     0.577350269189626d0, 0.577350269189626d0, 0.577350269189626d0  &
    ], [3,8])

    real(8), parameter :: np(3,8) = reshape([ &
    -1.0d0, -1.0d0,-1.0d0, &
     1.0d0, -1.0d0,-1.0d0, &
     1.0d0,  1.0d0,-1.0d0, &
    -1.0d0,  1.0d0,-1.0d0, &
    -1.0d0, -1.0d0, 1.0d0, &
     1.0d0, -1.0d0, 1.0d0, &
     1.0d0,  1.0d0, 1.0d0, &
    -1.0d0,  1.0d0, 1.0d0  &
    ], [3,8])

    integer(4), parameter :: C3D8_surf(4,6) = reshape([ &
        4, 3, 2, 1, &
        5, 6, 7, 8, &
        1, 2, 6, 5, &
        2, 3, 7, 6, &
        3, 4, 8, 7, &
        4, 1, 5, 8  ], [4,6])

    integer(4), parameter :: C3D8_edge(2,12) = reshape([ &
        1, 2, &
        2, 3, &
        3, 4, &
        4, 1, &
        5, 6, &
        6, 7, &
        7, 8, &
        8, 5, &
        1, 5, &
        2, 6, &
        3, 7, &
        4, 8  ], [2,12])

    ! public :: C3D8_integral_point
    ! public :: C3D8_get_global_deriv
    ! public :: C3D8_Bmat
    ! public :: C3D8_Dmat
    ! public :: C3D8_Kmat
    ! public :: get_interpolation_matrix_C3D8

contains

    subroutine C3D8_integral_point(i, r)
        implicit none
        integer(4), intent(in) :: i
        real(8), intent(out) :: r(3)
    
        r(1) = gsp(1,i)
        r(2) = gsp(2,i)
        r(3) = gsp(3,i)
    end subroutine C3D8_integral_point

    subroutine C3D8_get_global_deriv(node, r, dndx, det)
        implicit none
        real(8), intent(in) :: node(3,8)
        real(8), intent(in) :: r(3)
        real(8), intent(out) :: dndx(8,3)
        real(8), intent(out) :: det
        real(8) :: deriv(8,3), xj(3,3), inv(3,3)
        logical :: is_fail

        call C3D8_shapefunc_deriv(r, deriv)
        xj = matmul(node, deriv)
        call get_inverse_matrix_R_3d(xj, inv, det, is_fail)!//[ ] floating point
        dndx = matmul(deriv, inv)
    end subroutine C3D8_get_global_deriv

    subroutine C3D8_shapefunc(local, func)
        implicit none
        real(8), intent(in) :: local(3)
        real(8), intent(out) :: func(8)
    
        func(1) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0-local(3))
        func(2) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0-local(3))
        func(3) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0-local(3))
        func(4) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0-local(3))
        func(5) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0+local(3))
        func(6) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0+local(3))
        func(7) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0+local(3))
        func(8) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0+local(3))
    end subroutine C3D8_shapefunc

    subroutine C3D8_shapefunc_deriv(local, func)
        implicit none
        real(8), intent(in) :: local(3)
        real(8), intent(out) :: func(8,3)
    
        func(1,1) = -0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
        func(2,1) =  0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
        func(3,1) =  0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
        func(4,1) = -0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
        func(5,1) = -0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
        func(6,1) =  0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
        func(7,1) =  0.125d0*(1.0d0+local(2))*(1.0d0+local(3))
        func(8,1) = -0.125d0*(1.0d0+local(2))*(1.0d0+local(3))
    
        func(1,2) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
        func(2,2) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
        func(3,2) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
        func(4,2) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
        func(5,2) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(3))
        func(6,2) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
        func(7,2) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
        func(8,2) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(3))
    
        func(1,3) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
        func(2,3) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
        func(3,3) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
        func(4,3) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
        func(5,3) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
        func(6,3) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
        func(7,3) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
        func(8,3) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
    end subroutine C3D8_shapefunc_deriv

    subroutine get_inverse_matrix_R_3d(a, inv, det, is_fail)
        implicit none
        !> [in] 入力行列（サイズ [3, 3]）
        real(8), intent(in) :: a(3,3)
        !> [out] 逆行列（サイズ [3, 3]）
        real(8), intent(out) :: inv(3,3)
        !> [out] 行列式
        real(8), intent(out) :: det
        !> [out] 行列式が 0 以下であれば `true` となるフラグ
        logical, optional, intent(out) :: is_fail
        real(8) :: detinv
    
        if(present(is_fail)) is_fail = .false.
    
        det = a(1,1) * a(2,2) * a(3,3) &
            + a(2,1) * a(3,2) * a(1,3) &
            + a(3,1) * a(1,2) * a(2,3) &
            - a(3,1) * a(2,2) * a(1,3) &
            - a(2,1) * a(1,2) * a(3,3) &
            - a(1,1) * a(3,2) * a(2,3)
    
        if(det <= 0.0d0)then
          if(present(is_fail))then
            is_fail = .true.
          else
            call std_error_string("get_inverse_matrix_R_3d")
            call std_error_string("determinant is less than 0")
            call std_error_stop()
          endif
        endif
    
        detinv = 1.0d0/det
        inv(1,1) = detinv * ( a(2,2)*a(3,3) - a(3,2)*a(2,3))
        inv(1,2) = detinv * (-a(1,2)*a(3,3) + a(3,2)*a(1,3))
        inv(1,3) = detinv * ( a(1,2)*a(2,3) - a(2,2)*a(1,3))
        inv(2,1) = detinv * (-a(2,1)*a(3,3) + a(3,1)*a(2,3))
        inv(2,2) = detinv * ( a(1,1)*a(3,3) - a(3,1)*a(1,3))
        inv(2,3) = detinv * (-a(1,1)*a(2,3) + a(2,1)*a(1,3))
        inv(3,1) = detinv * ( a(2,1)*a(3,2) - a(3,1)*a(2,2))
        inv(3,2) = detinv * (-a(1,1)*a(3,2) + a(3,1)*a(1,2))
        inv(3,3) = detinv * ( a(1,1)*a(2,2) - a(2,1)*a(1,2))
    end subroutine get_inverse_matrix_R_3d

    subroutine get_inverse_matrix(n, a, inv)
        implicit none
        !> [in] 行列の大きさ
        integer(4), intent(in) :: n
        !> [in] 入力行列（サイズ [n, n]）
        real(8), intent(in) :: a(n,n)
        !> [out] 逆行列（サイズ [n, n]）
        real(8), intent(out) :: inv(n,n)
        integer(4) :: i, j, k
        real(8) :: b(n,n), tmp
    
        b = a
    
        inv = 0.0d0
        do i = 1, n
            inv(i,i) = 1.0d0
        enddo
    
        do i = 1, n
            if(b(i,i) == 0.0d0)then
                call std_error_string("get_inverse_matrix")
                call std_error_string("diagonal component is 0")
                call std_error_stop()
            endif
    
            tmp = 1.0d0/b(i,i)
    
            do j = 1, n
                b(j,i) =   b(j,i) * tmp
                inv(j,i) = inv(j,i) * tmp
            enddo
    
            do j = 1, n
                if(i /= j) then
                tmp = b(i,j)
                do k = 1, n
                    b(k,j) =   b(k,j) -   b(k,i) * tmp
                    inv(k,j) = inv(k,j) - inv(k,i) * tmp
                enddo
                endif
            enddo
        enddo
    end subroutine get_inverse_matrix

    subroutine C3D8_Bmat(dndx, u, B)
        implicit none
        integer(4) :: i, i1, i2, i3
        real(8) :: u(3,8), B(6,24), dndx(8,3), dudx(3,3)
    
        B = 0.0d0
        do i = 1,8
            i1 = 3*i-2
            i2 = 3*i-1
            i3 = 3*i
            B(1,i1) = dndx(i,1)
            B(2,i2) = dndx(i,2)
            B(3,i3) = dndx(i,3)
            B(4,i1) = dndx(i,2)
            B(4,i2) = dndx(i,1)
            B(5,i2) = dndx(i,3)
            B(5,i3) = dndx(i,2)
            B(6,i1) = dndx(i,3)
            B(6,i3) = dndx(i,1)
        enddo
    
        if(analysis_flag == 101)then
            dudx = 0.0d0
            dudx = matmul(u, dndx)
            do i = 1, 8
                i1 = 3*i-2
                i2 = 3*i-1
                i3 = 3*i
                B(1,i1) = B(1,i1) + dudx(1,1)*dndx(i,1)
                B(1,i2) = B(1,i2) + dudx(2,1)*dndx(i,1)
                B(1,i3) = B(1,i3) + dudx(3,1)*dndx(i,1)
                B(2,i1) = B(2,i1) + dudx(1,2)*dndx(i,2)
                B(2,i2) = B(2,i2) + dudx(2,2)*dndx(i,2)
                B(2,i3) = B(2,i3) + dudx(3,2)*dndx(i,2)
                B(3,i1) = B(3,i1) + dudx(1,3)*dndx(i,3)
                B(3,i2) = B(3,i2) + dudx(2,3)*dndx(i,3)
                B(3,i3) = B(3,i3) + dudx(3,3)*dndx(i,3)
                B(4,i1) = B(4,i1) + dudx(1,2)*dndx(i,1) + dudx(1,1)*dndx(i,2)
                B(4,i2) = B(4,i2) + dudx(2,2)*dndx(i,1) + dudx(2,1)*dndx(i,2)
                B(4,i3) = B(4,i3) + dudx(3,2)*dndx(i,1) + dudx(3,1)*dndx(i,2)
                B(5,i1) = B(5,i1) + dudx(1,2)*dndx(i,3) + dudx(1,3)*dndx(i,2)
                B(5,i2) = B(5,i2) + dudx(2,2)*dndx(i,3) + dudx(2,3)*dndx(i,2)
                B(5,i3) = B(5,i3) + dudx(3,2)*dndx(i,3) + dudx(3,3)*dndx(i,2)
                B(6,i1) = B(6,i1) + dudx(1,3)*dndx(i,1) + dudx(1,1)*dndx(i,3)
                B(6,i2) = B(6,i2) + dudx(2,3)*dndx(i,1) + dudx(2,1)*dndx(i,3)
                B(6,i3) = B(6,i3) + dudx(3,3)*dndx(i,1) + dudx(3,1)*dndx(i,3)
            enddo
        endif
    end subroutine C3D8_Bmat

    subroutine C3D8_Dmat(param, gauss, D)
        implicit none
        type(paramdef) :: param
        type(gaussdef) :: gauss
        real(8) :: D(6,6)

        call Dmat_elastic(param%E, param%nu, D)
    end subroutine C3D8_Dmat

    subroutine Dmat_elastic(E, mu, D)
        implicit none
        real(8) :: D(6,6), E, mu, g
    
        D = 0.0d0
        g = E / ((1.0d0+mu) * (1.0d0-2.0d0*mu))
    
        D(1,1) = g*(1.0d0-mu)
        D(1,2) = g*mu
        D(1,3) = g*mu
        D(2,1) = g*mu
        D(2,2) = g*(1.0d0-mu)
        D(2,3) = g*mu
        D(3,1) = g*mu
        D(3,2) = g*mu
        D(3,3) = g*(1.0d0-mu)
        D(4,4) = 0.5d0*g*(1.0d0-2.0d0*mu)
        D(5,5) = 0.5d0*g*(1.0d0-2.0d0*mu)
        D(6,6) = 0.5d0*g*(1.0d0-2.0d0*mu)
    end subroutine Dmat_elastic

    subroutine C3D8_update(mesh, var, param, icel, q)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        type(paramdef) :: param
        integer(4) :: i, in, icel
        real(8) :: x0(3,8), u(3,8), r(3), dndx(8,3), D(6,6), B(6,24)
        real(8) :: strain(6), stress(6), q(24), det
    
        q = 0.0d0
    
        do i = 1, 8
            in = mesh%elem(i,icel)
            x0(1,i) = mesh%node(1,in)
            x0(2,i) = mesh%node(2,in)
            x0(3,i) = mesh%node(3,in)
            u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
            u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
            u(3,i)  = var%u(3*in  ) + var%du(3*in  )
        enddo
    
        do i = 1, 8
            call C3D8_integral_point(i, r)
            call C3D8_get_global_deriv(x0, r, dndx, det)
            call C3D8_Bmat(dndx, u, B)
            call C3D8_get_starian(u, dndx, strain)
            var%gauss(i,icel)%strain = strain
        
            call Dmat_elastic(param%E, param%nu, D)
            var%gauss(i,icel)%stress = matmul(D, var%gauss(i,icel)%strain)

            q = q + matmul(var%gauss(i,icel)%stress, B)*det
        enddo
    end subroutine C3D8_update

    subroutine C3D8_get_starian(u, dndx, strain)
        integer(4) :: i, in, icel
        real(8) :: u(3,8), dndx(8,3), xj(3,3)
        real(8) :: strain(6)
    
        xj = matmul(u, dndx)
    
        strain(1) = xj(1,1)
        strain(2) = xj(2,2)
        strain(3) = xj(3,3)
        strain(4) =(xj(1,2) + xj(2,1))
        strain(5) =(xj(2,3) + xj(3,2))
        strain(6) =(xj(3,1) + xj(1,3))

        strain(1) = strain(1) + 0.5d0*dot_product(xj(:, 1), xj(:, 1))
        strain(2) = strain(2) + 0.5d0*dot_product(xj(:, 2), xj(:, 2))
        strain(3) = strain(3) + 0.5d0*dot_product(xj(:, 3), xj(:, 3))
        strain(4) = strain(4) + (xj(1,1)*xj(1,2) + xj(2,1)*xj(2,2) + xj(3,1)*xj(3,2))
        strain(5) = strain(5) + (xj(1,2)*xj(1,3) + xj(2,2)*xj(2,3) + xj(3,2)*xj(3,3))
        strain(6) = strain(6) + (xj(1,1)*xj(1,3) + xj(2,1)*xj(2,3) + xj(3,1)*xj(3,3))

    end subroutine C3D8_get_starian

    subroutine C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)
        implicit none
        type(vardef) :: var
        integer(4) :: i, j, k, icel
        real(8) :: inv(8,8)
        real(8) :: nstrain(8,6), nstress(8,6)
        real(8) :: estrain(6), estress(6)
    
        nstrain  = 0.0d0
        nstress  = 0.0d0
        estrain  = 0.0d0
        estress  = 0.0d0
    
        do i = 1, 8
            do j = 1, 8
                do k = 1, 6
                nstrain(i,k) = nstrain(i,k) + inv(i,j) * var%gauss(j,icel)%strain(k)
                nstress(i,k) = nstress(i,k) + inv(i,j) * var%gauss(j,icel)%stress(k)
                enddo
            enddo
        enddo
    
        do i = 1, 8
            do j = 1, 6
                estrain(j) = estrain(j) + var%gauss(i,icel)%strain(j)
                estress(j) = estress(j) + var%gauss(i,icel)%stress(j)
            enddo
        enddo
        estrain = estrain/8.0d0
        estress = estress/8.0d0
    end subroutine C3D8_get_nodal_values

    subroutine C3D8_Kmat(D, B, wg, det, stress, dndx, stiff)
        implicit none
        integer(4) :: i, j, k
        real(8) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det
        real(8) :: stress(6), S(9,9), BN(9,24), SBN(9,24), dndx(8,3)
    
        DB = matmul(D, B)
        do i = 1, 24
            do j = 1, 24
                do k = 1, 6
                    stiff(j,i) = stiff(j,i) + B(k,j)*DB(k,i)*wg*det
                enddo
            enddo
        enddo

        if(analysis_flag == 101)then
            BN = 0.0d0
            do j = 1, 8
                BN(1, 3*j-2) = dndx(j, 1)
                BN(2, 3*j-1) = dndx(j, 1)
                BN(3, 3*j  ) = dndx(j, 1)
                BN(4, 3*j-2) = dndx(j, 2)
                BN(5, 3*j-1) = dndx(j, 2)
                BN(6, 3*j  ) = dndx(j, 2)
                BN(7, 3*j-2) = dndx(j, 3)
                BN(8, 3*j-1) = dndx(j, 3)
                BN(9, 3*j  ) = dndx(j, 3)
            enddo

            S = 0.0d0
            do j = 1, 3
                S(j  , j  ) = stress(1)
                S(j  , j+3) = stress(4)
                S(j  , j+6) = stress(6)
                S(j+3, j  ) = stress(4)
                S(j+3, j+3) = stress(2)
                S(j+3, j+6) = stress(5)
                S(j+6, j  ) = stress(6)
                S(j+6, j+3) = stress(5)
                S(j+6, j+6) = stress(3)
            enddo

            SBN = matmul(S, BN)
            do i = 1, 24
                do j = 1, 24
                    do k = 1, 9
                        stiff(j,i) = stiff(j,i) + BN(k,j)*SBN(k,i)*wg*det
                    enddo
                enddo
            enddo
        endif
    end subroutine C3D8_Kmat
end module mod_element