module mod_analysis
    use mod_util
    use mod_io
    use mod_matrix
    use mod_solver
    use mod_update
    implicit none
    
contains
    subroutine small_deformation_analysis(mesh, param, var)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        type(vardef) :: var
        integer(4) :: i,j
        real(8),allocatable :: AT(:,:)
        character :: esc*1 = char(27)
        allocate(AT(ndof*mesh%nnode,ndof*mesh%nnode), source = 0.0d0)

        call initialize_mesh(mesh, var)
        call initialize_matrix(mesh)

        call load_condition(var, param)
        write(*,*)esc//"[32m"//'load condition is done.'//esc//"[0m"

        call get_stiff_matrix(mesh, param, var)
        write(*,*)esc//"[32m"//'generating matrix is done.'//esc//"[0m"
        call bound_condition(mesh, param, var)
        write(*,*)esc//"[32m"//'applying bound conditions is done.'//esc//"[0m"
        call get_RHS(mesh, var)

        call solver(mesh, var)
        write(*,*)esc//"[32m"//'solving linear equations is done.'//esc//"[0m"

        call delta_u_update(var)
        call stress_update(mesh, var, param)
        write(*,*)esc//"[32m"//'calculatig stress is done.'//esc//"[0m"

        call u_update(var)
        call outout_res(mesh, param, var)
        write(*,*)esc//"[32m"//'visualizing is done.'//esc//"[0m"
        call writeout_u(var, param, mesh)

    end subroutine small_deformation_analysis

    subroutine large_deformation_analysis(mesh, param, var)
        implicit none
        type(meshdef) :: mesh
        type(paramdef) :: param
        type(vardef) :: var
        integer(4) :: NRiter
        character :: esc*1 = char(27)

        call initialize_mesh(mesh, var)
        call initialize_matrix(mesh)

        call load_condition(var, param)

        write(*,*)esc//"[32m"//'load condition is done.'//esc//"[0m"
        write(*,*)esc//"[36m"//'<start NR loop>'//esc//"[0m"

        do NRiter = 1, param%max_nrstep
            write(*,'(i0,a)')NRiter,esc//"[36m"//'th NR loop'//esc//"[0m"
            call get_stiff_matrix(mesh, param, var)
            write(*,*)esc//"[32m"//'generating matrix is done.'//esc//"[0m"
            call get_RHS(mesh, var)
            call bound_condition(mesh, param, var)
            write(*,*)esc//"[32m"//'applying bc is done.'//esc//"[0m"

            if(is_convergence(mesh, var, NRiter)) then
                write(*,'(a,i0)')'The number of NR iteration is ',NRiter
                exit
            endif

            call solver(mesh, var)
            write(*,*)esc//"[32m"//'solving linear equations is done.'//esc//"[0m"

            call delta_u_update(var)
            call stress_update(mesh, var, param)
            write(*,*)esc//"[32m"//'calculatig stress is done.'//esc//"[0m"
        enddo

        call u_update(var)
        call outout_res(mesh, param, var)
        write(*,*)esc//"[32m"//'visualizing is done.'//esc//"[0m"
        call writeout_u(var, param, mesh)
    end subroutine large_deformation_analysis

end module mod_analysis