program main
    use mod_util
    use mod_io
    use mod_analysis

    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    type(vardef) :: var
    integer(4) :: ierr
    character :: esc*1 = char(27)

    call chdir('example', ierr)

    call input_param(param) !//[x] OK! input_param
    call input_mesh(mesh) !//[x] OK! input_mesh
    call input_bc(mesh, param) !//[x] OK! input_bc
    write(*,*)esc//"[32m"//'input is done.'//esc//"[0m"

    if(analysis_flag == 100) call small_deformation_analysis(mesh, param, var)
    if(analysis_flag == 101) call large_deformation_analysis(mesh, param, var)

    call chdir('../', ierr)
    write(*,*)esc//"[33m"//'all is done.'//esc//"[0m"

end program main