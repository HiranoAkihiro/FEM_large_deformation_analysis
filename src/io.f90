module mod_io
    use mod_util
    implicit none
    
contains

subroutine input_param(param)
    implicit none
    type(paramdef) :: param
    integer(4) :: i, n

    open(10,file="input.dat",status="old")
        read(10,*)param%E
        read(10,*)param%nu
        read(10,*)param%rho
        read(10,*)param%max_nrstep
        read(10,*)analysis_flag
    close(10)

end subroutine

subroutine input_mesh(mesh)
    implicit none
    type(meshdef) :: mesh
    integer(4) :: i, j

    open(11,file="node.dat",status="old")
        read(11,*)mesh%nnode
        allocate(mesh%node(ndof,mesh%nnode))
        do i = 1, mesh%nnode
            read(11,*)mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
        enddo
    close(11)

    open(12,file="elem.dat",status="old")
    read(12,*)mesh%nelem, mesh%nbase_func
        allocate(mesh%elem(mesh%nbase_func, mesh%nelem))
        do i = 1, mesh%nelem
            read(12,*)(mesh%elem(j,i), j = 1, mesh%nbase_func)
        enddo
    close(12)
end subroutine input_mesh

subroutine input_bc(mesh, param)
    implicit none
    type(meshdef) :: mesh
    type(paramdef) :: param
    integer(4) :: i, n_dof

    open(20, file = "bc.dat", status = "old")
        read(20,*) param%nbound, n_dof

        allocate(param%ibound(2, param%nbound))
        allocate(param%bound(param%nbound))

        do i = 1, param%nbound
            read(20,*) param%ibound(1,i), param%ibound(2,i), param%bound(i)
        enddo
    close(20)

    open(20, file = "load.dat", status = "old")
        read(20,*) param%ncload, n_dof

        allocate(param%icload(2, param%ncload))
        allocate(param%cload(param%ncload))

        do i = 1, param%ncload
            read(20,*) param%icload(1,i), param%icload(2,i), param%cload(i)
        enddo
    close(20)

end subroutine input_bc

    subroutine outout_res(mesh, param, var)
        implicit none
        type(paramdef) :: param
        type(meshdef) :: mesh
        type(vardef) :: var
        integer(4) :: i, id, nnode, nelem
        character :: cstep*5, cnum*5, output_dir*100

        !call soild_debug_header("outout_res")

        nnode = mesh%nnode
        nelem = mesh%nelem

        output_dir = "visual/"
        call system('if [ ! -d visual ]; then (echo "** create visual"; mkdir -p visual); fi')

        open(20, file='visual/u.dat', status='replace')
            write(20,"(i0)")nnode
            do i = 1, nnode
                write(20,"(1p3e22.14)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
            enddo
        close(20)

        call convert_to_real(mesh, var)

        write(cstep,"(i5.5)")0

        open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.pvtu', status='replace')
            write(20,"(a)")'<?xml version="1.0"?>'
            write(20,"(a)")'<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">'
            write(20,"(a)")'<PUnstructuredGrid>'
            write(20,"(a)")'<PPoints>'
            write(20,"(a)")'<PDataArray type="Float32" NumberOfComponents="3"/>'
            write(20,"(a)")'</PPoints>'
            write(20,"(a)")'<PCells>'
            write(20,"(a)")'<PDataArray type="Int32" Name="connectivity" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Int32" Name="offsets" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Int32" Name="types" format="appended"/>'
            write(20,"(a)")'</PCells>'
            write(20,"(a)")'<PPointData>'
            write(20,"(a)")'<PDataArray type="Float32" Name="disp" NumberOfComponents="3" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="nstress" NumberOfComponents="6" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="nmises" NumberOfComponents="1" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="appended"/>'
            write(20,"(a)")'</PPointData>'
            write(20,"(a)")'<PCellData>'
            write(20,"(a)")'<PDataArray type="Float32" Name="estrain" NumberOfComponents="6" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="estress" NumberOfComponents="6" format="appended"/>'
            write(20,"(a)")'<PDataArray type="Float32" Name="emises" NumberOfComponents="1" format="appended"/>'
            write(20,"(a)")'</PCellData>'
            do i = 0, 0
                write(cnum,"(i0)") i
                write(20,"(a)")'<Piece Source="./result.'//trim(cstep)//'.'//trim(cnum)//'.vtu"/>'
            enddo
            write(20,"(a)")'</PUnstructuredGrid>'
            write(20,"(a)")'</VTKFile>'
        close(20)


        write(cnum,"(i0)")0

        open(20, file=trim(output_dir)//'result.'//trim(cstep)//'.'//trim(cnum)//'.vtu', status='replace')
            write(20,"(a)")'<?xml version="1.0"?>'
            write(20,"(a)")'<VTKFile type="UnstructuredGrid" version="1.0">'
            write(20,"(a)")'<UnstructuredGrid>'
            write(20,"(a,i0,a,i0,a)")'<Piece NumberOfPoints="', nnode, '" NumberOfCells="', nelem, '">'
            write(20,"(a)")'<Points>'
            write(20,"(a)")'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
            do i = 1, nnode
                write(20,"(1p3e20.12)")mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
            enddo
  
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'</Points>'
            write(20,"(a)")'<Cells>'
            write(20,"(a)")'<DataArray type="Int32" Name="connectivity" format="ascii">'
            do i = 1, nelem
                write(20,"(8i8)")mesh%elem(1,i)-1, mesh%elem(2,i)-1, mesh%elem(3,i)-1, mesh%elem(4,i)-1, &
                                mesh%elem(5,i)-1, mesh%elem(6,i)-1, mesh%elem(7,i)-1, mesh%elem(8,i)-1
            enddo
  
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Int32" Name="offsets" format="ascii">'
            do i = 1, nelem
                write(20,"(x,i0,$)")8*i
            enddo
            write(20,*)""
  
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="UInt8" Name="types" format="ascii">'
            do i = 1, nelem
                write(20,"(i3,$)")12
            enddo
            write(20,*)""
  
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'</Cells>'
    
            write(20,"(a)")'<PointData>'
            write(20,"(a)")'<DataArray type="Float32" Name="disp" NumberOfComponents="3" format="ascii">'
            do i = 1, nnode
                write(20,"(1p3e12.4)")var%u(3*i-2), var%u(3*i-1), var%u(3*i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="nstrain" NumberOfComponents="6" format="ascii">'
            do i = 1, nnode
                write(20,"(1p6e12.4)")var%nstrain(1,i), var%nstrain(2,i), var%nstrain(3,i), &
                                    & var%nstrain(4,i), var%nstrain(5,i), var%nstrain(6,i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="nstress" NumberOfComponents="6" format="ascii">'
            do i = 1, nnode
            write(20,"(1p6e12.4)")var%nstress(1,i), var%nstress(2,i), var%nstress(3,i), &
                                & var%nstress(4,i), var%nstress(5,i), var%nstress(6,i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="nmises" NumberOfComponents="1" format="ascii">'
            do i = 1, nnode
                write(20,"(1p6e12.4)")var%nmises(i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="nreaction" NumberOfComponents="3" format="ascii">'
            do i = 1, nnode
                write(20,"(1p6e12.4)")var%f_reaction(3*i-2), var%f_reaction(3*i-1), var%f_reaction(3*i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'</PointData>'
  
            write(20,"(a)")'<CellData>'
            write(20,"(a)")'<DataArray type="Float32" Name="estrain" NumberOfComponents="6" format="ascii">'
            do i = 1, nelem
                write(20,"(1p6e12.4)")var%estrain(1,i), var%estrain(2,i), var%estrain(3,i), &
                                    & var%estrain(4,i), var%estrain(5,i), var%estrain(6,i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="estress" NumberOfComponents="6" format="ascii">'
            do i = 1, nelem
                write(20,"(1p6e12.4)")var%estress(1,i), var%estress(2,i), var%estress(3,i), &
                                    & var%estress(4,i), var%estress(5,i), var%estress(6,i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'<DataArray type="Float32" Name="emises" NumberOfComponents="1" format="ascii">'
            do i = 1, nelem
                write(20,"(1p6e12.4)")var%emises(i)
            enddo
            write(20,"(a)")'</DataArray>'
            write(20,"(a)")'</CellData>'
  
            write(20,"(a)")'</Piece>'
            write(20,"(a)")'</UnstructuredGrid>'
            write(20,"(a)")'</VTKFile>'
        close(20)
  
    end subroutine outout_res

    subroutine convert_to_real(mesh, var)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        integer(4) :: i, j, nnode, nelem
        real(8) :: thr
    
        nnode = mesh%nnode
        nelem = mesh%nelem
        thr = 1.0d-30
    
        do i = 1, nnode
            if(var%nmises(i) < thr) var%nmises(i) = 0.0d0
        enddo
    
        do i = 1, nnode
            do j = 1, 3
                if(dabs(var%u    (3*i-3+j)) < thr) var%u    (3*i-3+j) = 0.0d0
                if(dabs(var%f_reaction(3*i-3+j)) < thr) var%f_reaction(3*i-3+j) = 0.0d0
            enddo
        enddo
    
        do i = 1, nnode
            do j = 1, 6
                if(dabs(var%nstrain(j,i)) < thr) var%nstrain(j,i) = 0.0d0
                if(dabs(var%nstress(j,i)) < thr) var%nstress(j,i) = 0.0d0
            enddo
        enddo
    
        do i = 1, nelem
            if(dabs(var%emises(i)) < thr) var%emises(i) = 0.0d0
        enddo
    
        do i = 1, nelem
            do j = 1, 6
                if(dabs(var%estrain(j,i)) < thr) var%estrain(j,i) = 0.0d0
                if(dabs(var%estress(j,i)) < thr) var%estress(j,i) = 0.0d0
            enddo
        enddo
    end subroutine convert_to_real
end module mod_io