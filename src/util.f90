module mod_util
    use mod_monolis
    implicit none
    integer(4), parameter :: ndof = 3
    integer(4) :: analysis_flag
    
    type gaussdef
        integer(4) :: is_yield
        real(8) :: strain(6)
        real(8) :: stress(6)
        real(8) :: eq_pstrain
        real(8) :: eq_pstrain_back
        real(8) :: eq_pstrain_trial
    end type gaussdef

    type meshdef
        integer(4) :: nnode                     !> 節点数
        integer(4) :: nelem                     !> 要素数
        integer(4) :: nbase_func                !> 要素を構成する節点数（形状関数の数）
        real(8), allocatable :: node(:,:)       !> 節点座標
        integer(4), allocatable :: elem(:,:)    !> 要素コネクティビティ
    end type meshdef

    type paramdef
        !> for NR loop
        integer(4) :: cur_nrstep
        integer(4) :: max_nrstep
        integer(4) :: analysis_flag
        !> for boundary condition
        integer(4) :: nbound
        integer(4), allocatable :: ibound(:,:)
        real(8), allocatable :: bound(:)

        integer(4) :: ncload
        integer(4), allocatable :: icload(:,:)
        real(8), allocatable :: cload(:)

        !> for elast-plactis
        real(8), allocatable :: strain_table(:)
        real(8), allocatable :: stress_table(:)

        !> for material property
        real(8) :: E, nu, rho

        !>for iterative methods and preconditioning
        integer(4) :: iter
        integer(4) :: prec
    end type paramdef
    
    type vardef
        !> for analysis
        real(8), allocatable :: x(:)  !> solution vector of Ax = b
        real(8), allocatable :: b(:)  !> solution vector of Ax = b
        real(8), allocatable :: u(:)  !> displacement
        real(8), allocatable :: du(:) !> delta displacement
        real(8), allocatable :: q(:)  !> internal force
        real(8), allocatable :: f(:)  !> external force
        real(8), allocatable :: f_reaction(:) !> reaction force

        !> for results
        type(gaussdef), allocatable :: gauss(:,:)
        !> Nodal components
        real(8), allocatable :: nstrain(:,:)
        real(8), allocatable :: nstress(:,:)
        real(8), allocatable :: nmises(:)
        !> Elemental components
        real(8), allocatable :: estrain(:,:)
        real(8), allocatable :: estress(:,:)
        real(8), allocatable :: emises(:)
    end type vardef

    type dense_struct
        real(8), pointer :: A(:,:) => null()
        real(8), pointer :: X(:) => null()
        real(8), pointer :: B(:) => null()
        logical, allocatable :: nonzero_patern(:,:)
    end type dense_struct

    type sparse_struct
        real(8), pointer :: A(:) => null()
        real(8), pointer :: X(:) => null()
        real(8), pointer :: B(:) => null()
    end type sparse_struct

    type matdef
        integer(4) :: Iarray(100) = 0           !> 整数パラメータ
        real(8) :: Rarray(100) = 0.0d0          !> 実数パラメータ
        integer(4) :: N                         !> 内部自由度数
        integer(4) :: NP                        !> 全自由度数
        integer(4) :: NDOF                      !> 1 ブロックの自由度
        type(dense_struct)  :: dence            !> 密行列形式
        type(sparse_struct) :: sparse           !> 疎行列形式
        ! type(matdef_separated_CSR) :: SCSR  !> 行列構造体（セパレート CSR 構造）
        ! type(matdef_CSR) :: CSR             !> 行列構造体（CSR 構造）
    end type matdef

    type(monolis_structure) :: mat
    type(monolis_com) :: com

contains

    subroutine initialize_matrix(mesh)
        implicit none
        type(meshdef) :: mesh

        call monolis_get_nonzero_pattern_by_simple_mesh_R(mat, mesh%nnode, 8, ndof, mesh%nelem, mesh%elem)
    end subroutine

    subroutine initialize_mesh(mesh, var)
        implicit none
        type(meshdef) :: mesh
        type(vardef) :: var
        integer(4) :: i, j
    
        allocate(var%gauss(8, mesh%nelem))
        allocate(var%nstrain(6, mesh%nnode), source = 0.0d0)
        allocate(var%nstress(6, mesh%nnode), source = 0.0d0)
        allocate(var%nmises (mesh%nnode), source = 0.0d0)
        allocate(var%estrain(6, mesh%nelem), source = 0.0d0)
        allocate(var%estress(6, mesh%nelem), source = 0.0d0)
        allocate(var%emises (mesh%nelem), source = 0.0d0)
    
        allocate(var%u (ndof*mesh%nnode), source = 0.0d0)
        allocate(var%du(ndof*mesh%nnode), source = 0.0d0)
        allocate(var%q (ndof*mesh%nnode), source = 0.0d0)
        allocate(var%f (ndof*mesh%nnode), source = 0.0d0)
        allocate(var%f_reaction (ndof*mesh%nnode), source = 0.0d0)
        allocate(var%x (ndof*mesh%nnode), source = 0.0d0)
        allocate(var%b (ndof*mesh%nnode), source = 0.0d0)

    
        do i = 1, mesh%nelem
            do j = 1, 8
                var%gauss(j,i)%is_yield = 0
                var%gauss(j,i)%strain = 0.0d0
                var%gauss(j,i)%stress = 0.0d0
                var%gauss(j,i)%eq_pstrain = 0.0d0
                var%gauss(j,i)%eq_pstrain_back = 0.0d0
                var%gauss(j,i)%eq_pstrain_trial = 0.0d0
            enddo
        enddo
    end subroutine initialize_mesh

    subroutine clear_mat_value(mat, mesh)
        implicit none
        type(matdef), intent(inout) :: mat
        type(meshdef), intent(in) :: mesh

        allocate(mat%dence%A(ndof*mesh%nnode,ndof*mesh%nnode))
        allocate(mat%dence%B(ndof*mesh%nnode))
        allocate(mat%dence%X(ndof*mesh%nnode))
        allocate(mat%dence%nonzero_patern(ndof*mesh%nnode,ndof*mesh%nnode), source = .false.)

        mat%dence%A = 0.0d0
        mat%dence%B = 0.0d0
        mat%dence%X = 0.0d0

    end subroutine clear_mat_value

end module mod_util