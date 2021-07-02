subroutine test()
  use castem_hho
  use iso_fortran_env
  implicit none
  type(ElementGeometry) :: eg
  type(ElementDescription) :: ed
  character(len=255) :: l
!  double precision, dimension (:) :: ewk
  double precision, dimension (:), allocatable :: ewk
!  double precision, dimension (1) :: swk
  double precision, dimension (50) :: gradients
  double precision, dimension (6) :: gauss_pts
  double precision, dimension (3) :: gauss_wts
  double precision, dimension (200) :: ufaces
  double precision, dimension (200) :: faces_corr
  double precision, dimension (200) :: stress
  double precision, dimension (700) :: ktan
  double precision, dimension (700) :: res
  double precision, dimension (700) :: system
  double precision, dimension (700) :: internal_forces
  type(ExitStatus) :: r
  eg % dim_eucli = 2
  eg % num_vertices = 3
  eg % num_faces = 3
  allocate ( eg % num_vertices_per_face (3) )
  allocate ( eg % connectivity (6) )
  allocate ( eg % vertices_coordinates (6) )
  eg % num_vertices_per_face = (/ 2, 2, 2 /)
  eg % connectivity = (/ 0, 1, &
       2 , 1 , &
       2 , 0 /)
  eg % vertices_coordinates = (/ &
       ! x- coordinate
       0.0, 1.0, 0.0, &
       ! y- coordinate
       0.0, 0.0, 1.0 /)
  r = print_element_description(eg)
  call get_command_argument(1, l)
  r = get_element_description(ed, TRIM(l), &
       "tri_1_hdg_equ_ss_get_element_description", eg)

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
  r = compute_system(ed, ewk, system, ktan, res, eg)
  ! DECONDENSATE
  ! IN : faces correction
  r = decondensate(ed, ewk, faces_corr, eg)
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

