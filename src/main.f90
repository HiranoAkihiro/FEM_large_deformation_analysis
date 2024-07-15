program main
    use mod_util
    use mod_io
    use mod_analysis
    use mod_monolis
    use mod_debug

    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(4) :: ierr
    character :: esc*1 = char(27)
    real(kdouble) :: t1, t2

    call chdir('example', ierr)

    call monolis_global_initialize()
    call monolis_initialize(mat)
    call monolis_com_initialize_by_parted_files(com, monolis_mpi_get_global_comm(), &
        & MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, "node.dat")
    
    t1 = monolis_get_time()

    call input_param(param) !//[x] OK! input_param
    call input_mesh(mesh) !//[x] OK! input_mesh
    call input_bc(mesh, param) !//[x] OK! input_bc
    write(*,*)esc//"[32m"//'input is done.'//esc//"[0m"

    if(analysis_flag == 100) call small_deformation_analysis(mesh, param, var)
    if(analysis_flag == 101) call large_deformation_analysis(mesh, param, var)

    call chdir('../', ierr)
    write(*,*)esc//"[33m"//'all is done.'//esc//"[0m"

    t2 = monolis_get_time()
    call plot_time("total ", t2 - t1)

    call monolis_finalize(mat)
    call monolis_global_finalize()

end program main