subroutine test()
    use castem_hho
    use iso_fortran_env
    implicit none
    type(ElementGeometry) :: eg
    type(ElementDescription) :: ed

    ! ORDER OF UNKNOWS : EXEMPLE ORDER 2
    ! CELL :
    ! CELL : UCELL_ORDER_0_X -> 1
    ! CELL : UCELL_ORDER_1_X -> X
    ! CELL : UCELL_ORDER_2_X -> Y
    ! CELL : UCELL_ORDER_3_X -> X²
    ! CELL : UCELL_ORDER_4_X -> XY
    ! CELL : UCELL_ORDER_5_X -> Y²
    ! CELL : UCELL_ORDER_0_Y
    ! CELL : UCELL_ORDER_1_Y
    ! CELL : UCELL_ORDER_2_Y
    ! CELL : UCELL_ORDER_3_X
    ! CELL : UCELL_ORDER_4_Y
    ! CELL : UCELL_ORDER_5_X
    ! FACE 0 :
    ! F0 : UF0_ORDER_0_X
    ! F0 : UF0_ORDER_1_X
    ! F0 : UF0_ORDER_2_X
    ! F0 : UF0_ORDER_0_Y
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

    character(len=255) :: l
    !  double precision, dimension (:) :: ewk
    double precision stabilization_parameter
    double precision E
    double precision nu
    integer(c_int64_t), parameter :: material_data_max = 20
    double precision, dimension (material_data_max) :: material_data
    double precision, dimension (:), allocatable :: ewk
    double precision, dimension (:), allocatable :: ewk_fixed
    double precision, dimension (:), allocatable :: ewk_dynamic
    integer(c_int64_t), parameter :: nf_max = 20
    integer(c_int64_t), parameter :: quad_size_max_triangle = 20
    integer(c_int64_t), parameter :: dim_field = 2
    integer(c_int64_t), parameter :: dim_eucli = 2
    integer(c_int64_t), parameter :: order_max = 3 + 1
    integer(c_int64_t), parameter :: gauss_pt_max_size = quad_size_max_triangle * nf_max
    integer(c_int64_t), parameter :: gradients_size_max = 9 * gauss_pt_max_size
    integer(c_int64_t), parameter :: stress_max_size = 6 * gauss_pt_max_size
    integer(c_int64_t), parameter :: ufaces_size_max = order_max * nf_max * dim_field
!    double precision, dimension (1000) :: ewk
    !  double precision, dimension (1) :: swk
    double precision, dimension (gradients_size_max) :: gradients
    double precision, dimension (gauss_pt_max_size * dim_eucli) :: gauss_pts
    double precision, dimension (gauss_pt_max_size) :: gauss_wts
    double precision, dimension (ufaces_size_max) :: ufaces
    double precision, dimension (ufaces_size_max) :: faces_correction
    double precision, dimension (stress_max_size) :: stress
    double precision, dimension (6 * 6 * gauss_pt_max_size) :: ktan
    double precision, dimension (ufaces_size_max) :: residual
    double precision, dimension (700) :: system
    double precision, dimension (700) :: internal_forces
    type(ExitStatus) :: r
    eg % dim_eucli = 2
    eg % num_vertices = 4
    eg % num_faces = 4
    allocate ( eg % num_vertices_per_face (4) )
    allocate ( eg % connectivity (8) )
    allocate ( eg % vertices_coordinates (8) )
    eg % num_vertices_per_face = (/ 2, 2, 2 /)
    eg % connectivity = (/ 0, 1, &
            2 , 1 , &
            2 , 3 , &
            3 , 0 /)
    eg % vertices_coordinates = (/ &
            ! x- coordinate
            0.0, 1.0, 1.0, 0.0, &
            ! y- coordinate
            0.0, 0.0, 1.0, 1.0 /)
    r = print_element_description(eg)
    call get_command_argument(1, l)
    r = get_element_description(ed, TRIM(l), "qua4_1_1_hho_ss", eg)
    !

    write(*,*) "workspace size: ", ed %element_workspace_size
    allocate (ewk (ed %element_workspace_size))

    r = initialize_workspace(ed, ewk, eg)
    write(*,*) "elem worksapce: ", ewk
    r = get_gauss_data(ed, gauss_pts, gauss_wts, eg)
    write(*,*) "gauss pt: ", gauss_pts
    write(*,*) "gauss wt: ", gauss_wts
    ! COMPUTE GRADIENTS
    ! IN : ufaces
    ! OUT : gradients
    r = compute_gradients(ed, ewk, gradients, ufaces, eg)
    ! COMPUTE INTERNAL FORCES
    ! IN : ufaces + stress
    ! OUT : internal forces
    r = compute_internal_forces(ed, ewk, internal_forces, stress, ufaces, eg)
    ! COMPUTE SYSTEM
    ! IN : Ktan + Res
    ! OUT : condensated system
    r = compute_system(ed, ewk, system, ktan, residual, eg)
    ! DECONDENSATE
    ! IN : faces correction
    r = decondensate(ed, ewk, faces_correction, eg)
    !   r = compute_gradients(ed, eg, ewk, swk, ufaces)
    ! CALL GET_ENVIRONMENT_VARIABLE( NAME="OS", LENGTH=long )
    ! IF (long > 0) THEN
    !   ALLOCATE( CHARACTER(LEN=long) :: path )
    !   CALL GET_ENVIRONMENT_VARIABLE( NAME="OS", VALUE=path )
    !   write(ioimp,*) "OS=", path
    !   DEALLOCATE( path )


    deallocate( eg % num_vertices_per_face)
    deallocate( eg % connectivity)
    deallocate( eg % vertices_coordinates)
    deallocate( ewk )
end subroutine test

program main
    call test()
end program main

