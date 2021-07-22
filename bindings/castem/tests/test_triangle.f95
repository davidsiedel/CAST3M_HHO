subroutine test()
    use castem_hho
    use iso_fortran_env
    implicit none
    !    ---------------------------------------------------------------------------------------------------------------
    !    ORDERING INFO
    !    ---------------------------------------------------------------------------------------------------------------
    ! ORDER OF UNKNOWS : EXEMPLE ORDER FACE = 2, ORDER CELL = 2, EUCLI DIM = 2
    ! CELL :
    ! CELL : UCELL_ORDER_0_X -> 1
    ! CELL : UCELL_ORDER_1_X -> X
    ! CELL : UCELL_ORDER_2_X -> Y
    ! CELL : UCELL_ORDER_3_X -> X²
    ! CELL : UCELL_ORDER_4_X -> XY
    ! CELL : UCELL_ORDER_5_X -> Y²
    ! CELL : UCELL_ORDER_0_Y ...
    ! CELL : UCELL_ORDER_1_Y
    ! CELL : UCELL_ORDER_2_Y
    ! CELL : UCELL_ORDER_3_X
    ! CELL : UCELL_ORDER_4_Y
    ! CELL : UCELL_ORDER_5_X
    ! FACE 0 :
    ! F0 : UF0_ORDER_0_X -> 1
    ! F0 : UF0_ORDER_1_X -> S avec S le parametre S(X,Y)
    ! F0 : UF0_ORDER_2_X -> S²
    ! F0 : UF0_ORDER_0_Y ...
    ! F0 : UF0_ORDER_1_Y
    ! F0 : UF0_ORDER_2_Y
    ! FACE 1 :
    ! F1 : UF0_ORDER_0_X
    ! F1 : UF1_ORDER_1_X
    ! F1 : UF2_ORDER_2_X
    ! F1 : UF0_ORDER_0_Y
    ! F1 : UF1_ORDER_1_Y
    ! F1 : UF2_ORDER_2_Y
    ! FACE 2 :
    ! F2 : UF0_ORDER_0_X
    ! F2 : UF1_ORDER_1_X
    ! F2 : UF2_ORDER_2_X
    ! F2 : UF0_ORDER_0_Y
    ! F2 : UF1_ORDER_1_Y
    ! F2 : UF2_ORDER_2_Y
    ! ...
    !    ---------------------------------------------------------------------------------------------------------------
    !    GENERIC MAXIMAL DATA IN 2D
    !    ---------------------------------------------------------------------------------------------------------------
    integer(c_int64_t), parameter :: nf_max = 20
    integer(c_int64_t), parameter :: quad_size_max_triangle = 12
    integer(c_int64_t), parameter :: dim_field = 2
    integer(c_int64_t), parameter :: dim_eucli = 2
    integer(c_int64_t), parameter :: order_max = 2
    integer(c_int64_t), parameter :: cl_max = 6
    integer(c_int64_t), parameter :: cr_max = 6
    integer(c_int64_t), parameter :: fk_max = 3
    integer(c_int64_t), parameter :: face_block_max = fk_max * nf_max * dim_field ! =
    integer(c_int64_t), parameter :: cell_block_max = cl_max * dim_field ! =
    integer(c_int64_t), parameter :: elem_block_max = cell_block_max + face_block_max ! =
    integer(c_int64_t), parameter :: B_matrix_max_size = elem_block_max * 9 ! =
    integer(c_int64_t), parameter :: Z_max_size = elem_block_max * elem_block_max ! =
    !    ---------------------------------------------------------------------------------------------------------------
    !    INITIALIZE TEST DATA
    !    ---------------------------------------------------------------------------------------------------------------
    character(len=255) :: l
    double precision stabilization_parameter
    double precision E
    double precision nu
    double precision, dimension (B_matrix_max_size) :: B_data
    double precision, dimension (Z_max_size) :: Z_data
    double precision, dimension (3) :: gauss_point_data
    double precision, dimension (1) :: gauss_weight_data
    integer(c_int64_t), parameter :: index = 1
    integer(c_int64_t) :: iloop = 1
    double precision, dimension (4) :: mat = (/ 0.5, 0.0, 0.0, 0.5/)
    integer(c_int64_t), parameter :: mat_dim = 2
    !    ---------------------------------------------------------------------------------------------------------------
    !    INITIALIZE OBJECTS
    !    ---------------------------------------------------------------------------------------------------------------
    type(ExitStatus) :: r
    type(ElementGeometry) :: eg
    type(ElementDescription) :: ed
    type(ElementFunctions) :: ef
    type(GenericFunctions) :: gf
    !    ---------------------------------------------------------------------------------------------------------------
    !    INITIALIZE GEOMETRY
    !    ---------------------------------------------------------------------------------------------------------------
    integer(c_int64_t), dimension(40) :: conn = (/      0, 1, 1, 2, 2, 0, 0, 0, 0, 0, &
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer(c_int64_t), dimension(40) :: verts = (/     0., 1., 0., 0., 0., 1., 0., 0., 0., 0., &
                                                        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
                                                        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., &
                                                        0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
    eg % connectivity = conn
    eg % vertices_coordinates = verts
    !    ---------------------------------------------------------------------------------------------------------------
    !    FOR ELEMENT LOOP
    !    ---------------------------------------------------------------------------------------------------------------
    ! call the library and build fotran objects
    call get_command_argument(1, l)
    r = get_element_description(ed, TRIM(l), "hho_1_1_sp_tri3_get_element_description")
    r = get_element_functions(ef, TRIM(l), "hho_1_1_sp_tri3_get_element_functions")
    r = get_generic_functions(gf, TRIM(l), "hho_get_generic_functions")
    ! print element geometry attributes
    write(*,*) "connectivity : ", eg % connectivity
    write(*,*) "vertices_coordinates : ", eg % vertices_coordinates
    ! print element description attributes
    write(*,*) "dim_eucli : ", ed % dim_eucli
    write(*,*) "num_vertices : ", ed % num_vertices
    write(*,*) "num_faces : ", ed % num_faces
    write(*,*) "num_vertices_per_face : ", ed % num_vertices_per_face
    write(*,*) "num_quadrature_points : ", ed % num_quadrature_points
    write(*,*) "dir_dof_face_unknown : ", ed % dir_dof_face_unknown
    write(*,*) "dir_dof_cell_unknown : ", ed % dir_dof_cell_unknown
    write(*,*) "dir_dof_gradient : ", ed % dir_dof_gradient
    write(*,*) "dim_field : ", ed % dim_field
    write(*,*) "dim_space_element : ", ed % dim_space_element
    write(*,*) "dim_face_block : ", ed % dim_face_block
    write(*,*) "dim_cell_block : ", ed % dim_cell_block
    write(*,*) "dim_MB : ", ed % dim_MB
    write(*,*) "dim_MB_matrices : ", ed % dim_MB_matrices
    write(*,*) "dim_MSTAB : ", ed % dim_MSTAB
    write(*,*) "dim_MKCC : ", ed % dim_MKCC
    write(*,*) "dim_MKCF : ", ed % dim_MKCF
    write(*,*) "dim_MVC : ", ed % dim_MVC
    ! get and print stabilization operator
    r = get_stabilization_operator(ef, eg, Z_data)
    write(*,*) "Z_data : ", Z_data ! sorted row-wise
    ! loop over indices of quadrature points
    do iloop = 1, ed % num_quadrature_points
        write(*,*) "iloop : ", iloop
        r = get_gradient_operator(ef, eg, B_data, iloop) ! sorted row-wise
        write(*,*) "B_data : ", B_data
        r = get_gauss_weight(ef, eg, gauss_weight_data, iloop)
        write(*,*) "gauss_weight_data : ", gauss_weight_data
        r = get_gauss_point(ef, eg, gauss_point_data, iloop)
        write(*,*) "gauss_point_data : ", gauss_point_data
    end do
    ! invert matrix
    r = invert_matrix(gf, mat, mat_dim)
    write(*,*) "mat : ", mat

end subroutine test

program main
    call test()
end program main

