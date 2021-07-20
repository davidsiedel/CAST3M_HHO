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
    !    MATERIAL DATA
    !    ---------------------------------------------------------------------------------------------------------------
    character(len=255) :: l
    double precision stabilization_parameter
    double precision E
    double precision nu
    !    ---------------------------------------------------------------------------------------------------------------
    !    GENERIC MAXIMAL DATA IN 2D
    !    ---------------------------------------------------------------------------------------------------------------
    integer(c_int64_t), parameter :: nf_max = 20
    integer(c_int64_t), parameter :: quad_size_max_triangle = 12
    integer(c_int64_t), parameter :: dim_field = 2
    integer(c_int64_t), parameter :: dim_eucli = 2
    integer(c_int64_t), parameter :: order_max = 2
    integer(c_int64_t), parameter :: elem_block_max = 6
    integer(c_int64_t), parameter :: cl_max = 6
    integer(c_int64_t), parameter :: cr_max = 6
    integer(c_int64_t), parameter :: fk_max = 3
    integer(c_int64_t), parameter :: face_block_max = fk_max * nf_max * dim_field ! =
    integer(c_int64_t), parameter :: cell_block_max = cl_max * dim_field ! =
    integer(c_int64_t), parameter :: elem_block_max = cell_block_max + face_block_max ! =
    integer(c_int64_t), parameter :: gauss_pt_max_size = quad_size_max_triangle * nf_max ! =
    integer(c_int64_t), parameter :: gradients_size_max = 9 * gauss_pt_max_size ! =
    integer(c_int64_t), parameter :: stresses_max_size = 6 * gauss_pt_max_size ! =
    integer(c_int64_t), parameter :: tangent_ops_max_size = 6 * 6 * gauss_pt_max_size ! =
    integer(c_int64_t), parameter :: B_matrices_max_size = gauss_pt_max_size * elem_block_max * 6 ! =
    integer(c_int64_t), parameter :: Z_max_size = elem_block_max * elem_block_max ! =
    !    ---------------------------------------------------------------------------------------------------------------
    !    GAUSS DATA
    !    ---------------------------------------------------------------------------------------------------------------
    double precision, dimension (gauss_pt_max_size * dim_eucli) :: gauss_pts
    double precision, dimension (gauss_pt_max_size) :: gauss_wts
    !    ---------------------------------------------------------------------------------------------------------------
    !    UNKNOWNS DATA AND CORRECTION DATA (ON FACES ONLY, SINCE GLOBAL SYSTEM DEPENDS ON FACES ONLY)
    !    ---------------------------------------------------------------------------------------------------------------
    double precision, dimension (face_block_max) :: faces_unknowns
    double precision, dimension (cell_block_max) :: cell_unknowns
    double precision, dimension (face_block_max) :: faces_correction
    !    ---------------------------------------------------------------------------------------------------------------
    !    GAUSS POINT COMPUATIONAL DATA
    !    ---------------------------------------------------------------------------------------------------------------
    double precision, dimension (9) :: gradient
    double precision, dimension (6) :: stress
    double precision, dimension (36) :: tangent_op
    double precision, dimension (elem_block_max * elem_block_max) :: K_elem
    double precision, dimension (elem_block_max) :: R_elem
    double precision, dimension (elem_block_max) :: F_int
    double precision, dimension (elem_block_max) :: F_ext
    !    ---------------------------------------------------------------------------------------------------------------
    !    INITIALIZE OBJECTS
    !    ---------------------------------------------------------------------------------------------------------------
    type(ExitStatus) :: r
    type(ElementGeometry) :: eg
    type(ElementDescription) :: ed
    type(ElementFunctions) :: ef
    double precision, dimension (B_matrices_max_size) :: B_matrices
    double precision, dimension (Z_max_size) :: Z
    !    ---------------------------------------------------------------------------------------------------------------
    !    FOR ELEMENT LOOP
    !    ---------------------------------------------------------------------------------------------------------------
    do elem = 0, 10
        call get_command_argument(1, l)
        r = get_element_description(ed, TRIM(l), "qua4_1_1_hho_ss")
        r = get_element_functions(ef, TRIM(l), "qua4_1_1_hho_ss")
        !    -----------------------------------------------------------------------------------------------------------
        !    FILL ELEMENT GEOMETRY
        !    -----------------------------------------------------------------------------------------------------------
        !    !!!!!!!!!!!!!!!!!!! WRITE FULL ARRAYS ? (SINCE MAX SIZE IS 40)
        eg % connectivity = (/  0 , 1, &
                2 , 1 , &
                2 , 3 , &
                3 , 0 /)
        eg % vertices_coordinates = (/ &
                ! x- coordinate
                0.0, 1.0, 1.0, 0.0, &
                ! y- coordinate
                0.0, 0.0, 1.0, 1.0 /)
        r = compute_operators(ef, B_matrices, Z, eg)
        !    -----------------------------------------------------------------------------------------------------------
        !    COPY OPERATORS DATA TO CHAMELEM
        !    -----------------------------------------------------------------------------------------------------------
    end do
    !    ---------------------------------------------------------------------------------------------------------------
    !    FOR ELEMENT LOOP -> COMPUTATION STEP
    !    ---------------------------------------------------------------------------------------------------------------
    do elem = 0, 10
        call get_command_argument(1, l)
        r = get_element_description(ed, TRIM(l), "qua4_1_1_hho_ss")
        r = get_element_functions(ef, TRIM(l), "qua4_1_1_hho_ss")
        r = get_gauss_data(ef, gauss_pts, gauss_wts, eg)
        do index = 0, ed % num_quadrature_points
            r = compute_gradients(ef, B_matrices, Z, gradient, faces_unknowns, cell_unknowns, index)
        end do
        !    -----------------------------------------------------------------------------------------------------------
        !    INTREGRATE BEHAVIOUR LAW AND COMPUTE STRESS, TANGENT OPERATOR
        !    -----------------------------------------------------------------------------------------------------------
        do index = 0, ed % num_quadrature_points
            r = compute_internal_forces(ef, B_matrices, Z, F_int, stress, stabilization_parameter, faces_unknowns, cell_unknowns, index)
            r = compute_tangent_matrix(ef, B_matrices, Z, K_elem, tangent_op, stabilization_parameter, faces_unknowns, cell_unknowns, index)
        end do
        !    -----------------------------------------------------------------------------------------------------------
        !    COMPUTE EXTERNAL FORCES AND DEDUCE R_elem
        !    -----------------------------------------------------------------------------------------------------------
        r = condensate_system(ef, K_elem, R_elem)
        !    -----------------------------------------------------------------------------------------------------------
        !    ASSEMBLE AND SOLVE GLOBAL SYSTEM
        !    -----------------------------------------------------------------------------------------------------------
        r = decondensate(ef, ewk, faces_correction, eg)
    end do

end subroutine test

program main
    call test()
end program main

