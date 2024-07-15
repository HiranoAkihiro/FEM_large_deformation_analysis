module mod_update
    use mod_util
    use mod_element
    implicit none
    
contains

    subroutine delta_u_update(var)
        implicit none
        type(vardef) :: var

        !call soild_debug_header("delta_u_update")
        var%du = var%du + var%X
    end subroutine delta_u_update

    subroutine u_update(var)
        implicit none
        type(vardef) :: var
    
        !call soild_debug_header("u_update")
        var%u = var%u + var%du
    end subroutine u_update

    subroutine stress_update(mesh, var, param)
        use mod_matrix
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        type(paramdef) :: param
        integer(4) :: i, j, in, icel
        real(8) :: func(8,8), inv(8,8), tmp
        real(8) :: nstrain(8,6), nstress(8,6)
        real(8) :: estrain(6),   estress(6)
        real(8) :: q(24), r(3)
        integer(4), allocatable :: inode(:)
    
        !call soild_debug_header("stress_update")
    
        call init_nodal_strain_and_stress(mesh, var)
        call get_interpolation_matrix_C3D8(inv)
    
        allocate(inode(mesh%nnode), source = 0)
        var%q = 0.0d0
    
        do icel = 1, mesh%nelem
            call C3D8_update(mesh, var, param, icel, q)
            call C3D8_get_nodal_values(var, icel, inv, nstrain, nstress, estrain, estress)
        
            do i = 1, 8
                in = mesh%elem(i,icel)
                inode(in) = inode(in) + 1
                do j = 1, 6
                    var%nstrain(j,in) = var%nstrain(j,in) + nstrain(i,j)
                    var%nstress(j,in) = var%nstress(j,in) + nstress(i,j)
                enddo
                var%q(3*in-2) = var%q(3*in-2) + q(3*i-2)
                var%q(3*in-1) = var%q(3*in-1) + q(3*i-1)
                var%q(3*in  ) = var%q(3*in  ) + q(3*i  )
            enddo
    
            do j = 1, 6
                var%estrain(j,icel) = estrain(j)
                var%estress(j,icel) = estress(j)
            enddo
        enddo
    
        !> get average of nodal strain and stress
        do i = 1, mesh%nnode
            tmp = 1.0d0/dble(inode(i))
            do j = 1, 6
                var%nstrain(j,i) = var%nstrain(j,i) * tmp
                var%nstress(j,i) = var%nstress(j,i) * tmp
            enddo
        enddo
    
        !> get nodal Mises stress
        do i = 1, mesh%nnode
            call get_mises(var%nstress(1:6,i), var%nmises(i))
        enddo
    
        !> get elemental Mises stress
        do i = 1, mesh%nelem
            call get_mises(var%estress(1:6,i), var%emises(i))
        enddo
    end subroutine stress_update

    subroutine init_nodal_strain_and_stress(mesh, var)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        integer(4) :: i, j
        do i = 1, mesh%nnode
            do j = 1, 6
                var%nstrain(j,i) = 0.0d0
                var%nstress(j,i) = 0.0d0
            enddo
            enddo
    end subroutine init_nodal_strain_and_stress

    subroutine get_interpolation_matrix_C3D8(inv)
        implicit none
        integer(4) :: i
        real(8) :: func(8,8), inv(8,8), r(3)
    
        do i = 1, 8
            call C3D8_integral_point(i, r)
            call C3D8_shapefunc(r, func(i,:))
        enddo
        call get_inverse_matrix(8, func, inv)
    end subroutine get_interpolation_matrix_C3D8

    subroutine get_mises(s, mises)
        implicit none
        real(8) :: mises, s(6)
        real(8) :: s11, s22, s33, s12, s23, s13, ps, smises
    
        s11 = s(1)
        s22 = s(2)
        s33 = s(3)
        s12 = s(4)
        s23 = s(5)
        s13 = s(6)
        ps = (s11 + s22 + s33) / 3.0d0
        smises = 0.5d0 * ((s11-ps)**2 + (s22-ps)**2 + (s33-ps)**2) + s12**2 + s23**2 + s13**2
        mises  = dsqrt( 3.0d0 * smises )
    end subroutine get_mises
end module mod_update