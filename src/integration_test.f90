program integration_test
    use mod_util
    use mod_io
    use mod_analysis

    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(4) :: ierr
    character :: esc*1 = char(27)

    call chdir('directory_for_test', ierr)

    call input_param(param)
    call input_mesh(mesh)
    call input_bc(mesh, param)
    write(*,*)esc//"[32m"//'input is done.'//esc//"[0m"

    if(analysis_flag == 100) call small_deformation_analysis(mesh, param, var)
    if(analysis_flag == 101) call large_deformation_analysis(mesh, param, var)

    call chdir('../', ierr)
    write(*,*)esc//"[33m"//'program is done.'//esc//"[0m"

end program integration_test