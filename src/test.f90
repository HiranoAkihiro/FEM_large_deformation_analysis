program test
    use mod_util
    use mod_debug
    use mod_element
    use mod_update
    use mod_io
    use mod_matrix
    use mod_solver
    use mod_analysis
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(4) :: ierr
    integer(4) :: i, in, icel
    integer(4) :: elem(8)
    real(8) :: x(3,8), stiff(24,24), u(3,8)
    real(8) :: r(3), wg, det
    real(8) :: B(6,24), D(6,6), dndx(8,3)
    character :: esc*1 = char(27)

    call chdir('directory_for_test', ierr)

    call input_param(param)
    call input_mesh(mesh)
    call input_bc(mesh, param)
    write(*,*)esc//"[32m"//'input is done.'//esc//"[0m"
    !> lda
    call clear_mat_value(mat, mesh)
    call initialize_mesh(mesh, var)
    call load_condition(var, param)
    write(*,*)esc//"[32m"//'load condition is done.'//esc//"[0m"
    write(*,*)esc//"[36m"//'<start NR loop>'//esc//"[0m"

    do icel = 1, mesh%nelem
        call get_element_node_id(icel, mesh%elem, elem)
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
            call C3D8_integral_point(i, r)
            
        enddo
    enddo


end program test