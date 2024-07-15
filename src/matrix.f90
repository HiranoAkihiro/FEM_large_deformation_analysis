module mod_matrix
    use mod_util
    use mod_element
    implicit none
    
contains

    subroutine get_stiff_matrix(mesh, param, var)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        type(vardef) :: var
        integer(4) :: i, icel
        integer(4) :: elem(8)
        real(8) :: stiff(24,24), x(3,8)

        call monolis_clear_mat_value_R(mat)

        do icel = 1, mesh%nelem
            call get_element_node_id(icel, mesh%elem, elem)
            call C3D8_stiff(mesh, var, param, icel, stiff)
            call monolis_add_matrix_to_sparse_matrix_R(mat, 8, elem, stiff)
            ! call add_element_matrix_to_dense_matrix(mat, 8, elem, stiff) !//[ ] add_element_matrix_to_dense_matrix
        enddo
    end subroutine get_stiff_matrix

    subroutine C3D8_stiff(mesh, var, param, icel, stiff)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        type(paramdef) :: param
        integer(4) :: i, in, icel
        integer(4) :: elem(8)
        real(8) :: x(3,8), stiff(24,24), u(3,8)
        real(8) :: r(3), wg, det
        real(8) :: B(6,24), D(6,6), dndx(8,3)

        wg    = 1.0d0
        stiff = 0.0d0

        do i = 1, 8
            in = mesh%elem(i,icel)
            x(1,i) = mesh%node(1,in)
            x(2,i) = mesh%node(2,in)
            x(3,i) = mesh%node(3,in)
            u(1,i)  = var%u(3*in-2) + var%du(3*in-2)
            u(2,i)  = var%u(3*in-1) + var%du(3*in-1)
            u(3,i)  = var%u(3*in  ) + var%du(3*in  )
        enddo

        do i = 1, 8
            call C3D8_integral_point(i, r) ! //[x] OK! C3D8_integral_point
            call C3D8_get_global_deriv(x, r, dndx, det)! //[x] OK! C3D8_get_global_deriv
            call C3D8_Bmat(dndx, u, B) ! //[ ] monolis C3D8_Bmat
            call C3D8_Dmat(param, var%gauss(i,icel), D) !//[ ] monolis C3D8_Dmat
            call C3D8_Kmat(D, B, wg, det, var%gauss(i,icel)%stress, dndx, stiff) ! //[ ] monolis C3D8_Kmat
        enddo

    end subroutine C3D8_stiff

    subroutine load_condition(var, param)
        implicit none
        type(paramdef) :: param
        type(vardef) :: var
        integer(4) :: i, in, dof
        real(8) :: val

        var%f = 0.0d0
        do i = 1, param%ncload
            in = param%icload(1, i)
            dof = param%icload(2, i)
            val = param%cload(i)
            if(ndof < dof) stop "*** error: 3 < dof"
                var%f(ndof*(in-1) + dof) = val
        enddo
    end subroutine load_condition

    subroutine get_element_node_id(eid, elem, elemid)
        implicit none
        integer(4) :: i, eid, elem(:,:), elemid(:)
        do i = 1, 8
            elemid(i) = elem(i,eid)
        enddo
    end subroutine get_element_node_id

    subroutine add_element_matrix_to_dense_matrix(mat, n_base, connectivity, e_mat)
        implicit none
        type(matdef) :: mat
        integer(4), intent(in) :: n_base !> 要素を構成する節点数
        integer(4), intent(in) :: connectivity(n_base)
        integer(4) :: i, j, k, l
        real(8), intent(in) :: e_mat(:,:) !> 要素剛性行列

        do l=1, n_base
            do k=1, n_base
                do j=1, ndof
                    do i=1, ndof
                        mat%dence%A(ndof*(connectivity(k)-1)+i, ndof*(connectivity(l)-1)+j) &
                        = e_mat(ndof*(k-1)+i, ndof*(l-1)+j)
                        mat%dence%nonzero_patern(ndof*(connectivity(k)-1)+i, ndof*(connectivity(l)-1)+j) &
                        = .true.
                    enddo
                enddo
            enddo
        enddo
    end subroutine add_element_matrix_to_dense_matrix

    subroutine get_RHS(mesh, var)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var

        var%B = var%f - var%q
    end subroutine get_RHS

    subroutine bound_condition(mesh, param, var)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        type(vardef) :: var
        integer(4) :: i, in, dof, nb
        ! integer(4), allocatable :: indexR(:), itemR(:), permA(:)
        real(8) :: val

        do nb = 1, param%nbound
            in  = param%ibound(1, nb)
            dof = param%ibound(2, nb)
            val = param%bound(nb) - var%u(ndof*(in-1) + dof) - var%du(ndof*(in-1) + dof)
            if(ndof < dof) stop "*** error: 3 < dof"
            call monolis_set_Dirichlet_bc_R(mat, var%B, in, dof, val)
        enddo
    end subroutine bound_condition

    subroutine set_Dirichlet_bc(mat, B, node_id, ndof_bc, val)
        implicit none
        type(matdef), intent(inout) :: mat
        real(8), intent(inout) :: B(:)
        integer(4), intent(in) :: node_id
        integer(4), intent(in) :: ndof_bc
        integer(4) :: i
        real(8), intent(in) :: val

        mat%dence%A(:, ndof*(node_id-1)+ndof_bc) = 0.0d0
        mat%dence%A(ndof*(node_id-1)+ndof_bc, :) = 0.0d0
        mat%dence%A(ndof*(node_id-1)+ndof_bc, ndof*(node_id-1)+ndof_bc) = 1.0d0
        mat%dence%nonzero_patern(:, ndof*(node_id-1)+ndof_bc) = .false.
        mat%dence%nonzero_patern(ndof*(node_id-1)+ndof_bc, :) = .false.
        mat%dence%nonzero_patern(ndof*(node_id-1)+ndof_bc, ndof*(node_id-1)+ndof_bc) = .true.

        ! ここで考える片持ち梁ではディリクレ境界条件は強制変位０なので、右辺へのは操作は無しとした。

    end subroutine set_Dirichlet_bc

    subroutine convert_dense_to_sparse(mat)
        implicit none
        type(matdef) :: mat

        ! call convert_dense_to_sparse_main(mat%dence%A)

    end subroutine convert_dense_to_sparse
end module mod_matrix