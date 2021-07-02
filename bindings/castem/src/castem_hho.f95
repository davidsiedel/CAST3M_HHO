module castem_hho
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_ptr, c_funptr, c_char, c_null_ptr, c_loc
  use, intrinsic :: iso_fortran_env, only: int64
  use castem_hho_fortran_utilities
  !
  type :: ElementGeometry
     integer(c_int64_t) :: dim_eucli
     integer(c_int64_t) :: num_vertices
     integer(c_int64_t) :: num_faces
     integer(c_int64_t), dimension (:), pointer :: num_vertices_per_face
     integer(c_int64_t), dimension (:), pointer :: connectivity
     double precision, dimension (:), pointer :: vertices_coordinates
  end type ElementGeometry
  !
  type, bind(c) :: ElementDescription
     integer(c_int64_t) element_workspace_size
     integer(c_int64_t) num_cell_quadrature_points
     integer(c_int64_t) element_size
     integer(c_int64_t) face_size
     type(c_funptr) get_gauss_data
     type(c_funptr) initialize_workspace
     type(c_funptr) compute_gradients
     type(c_funptr) compute_internal_forces
     type(c_funptr) compute_system
  end type ElementDescription
  !
  type, bind(c) :: CElementGeometry
     integer(c_int64_t) :: dim_eucli
     integer(c_int64_t) :: num_vertices
     integer(c_int64_t) :: num_faces
     type(c_ptr) :: num_vertices_per_face
     type(c_ptr) :: connectivity
     type(c_ptr) :: vertices_coordinates
  end type CElementGeometry
  ! enumeration
  enum , bind(c) 
     enumerator :: CASTEM_HHO_SUCCESS = 0
     enumerator :: CASTEM_HHO_FAILURE = 1
  end enum
  type, bind(c) :: ExitStatus
     integer(c_int) :: ExitStatus
     type(c_ptr) :: msg = c_null_ptr
  end type ExitStatus
  interface
     function report_success() bind(c,name = 'castem_hho_report_success') result(r)
       import ExitStatus
       implicit none
       type(ExitStatus) :: r
     end function report_success
  end interface
  !
contains
  !
  function convert_element_geometry(d) result(r)
    implicit none
    type(ElementGeometry), intent(in) :: d
    type(CElementGeometry) :: r
    r % dim_eucli = d % dim_eucli
    r % num_vertices = d % num_vertices
    r % num_faces = d % num_faces
    r % num_vertices_per_face = c_loc(d % num_vertices_per_face)
    r % connectivity = c_loc(d % connectivity)
    r % vertices_coordinates = c_loc(d % vertices_coordinates)
  end function convert_element_geometry
  !
  function report_failure(msg) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_loc
    implicit none
    interface
       function report_failure_wrapper(msg) bind(c,name = 'castem_hho_report_failure') result(r)
         import ExitStatus, c_ptr, c_char
         implicit none
         character(len=1,kind=c_char), dimension(*), intent(in) :: msg
         type(ExitStatus) :: r
       end function report_failure_wrapper
    end interface
    character(len=*), intent(in) :: msg
    type(ExitStatus) :: r
    character(len=1,kind=c_char) :: s(len(msg)+1)
    s = convert_fortran_string(msg)
    r = report_failure_wrapper(s)
  end function report_failure
  !
  function get_error_message(s) result(msg)
    implicit none
    type(ExitStatus) :: s
    character(len=:), allocatable :: msg
    msg = convert_c_string(s%msg)
  end function get_error_message
  !
  function print_element_description(gd) result(r)
    implicit none
    interface
       function print_element_description_wrapper(gd) &
            bind(c,name = 'castem_hho_print_element_geometry') &
            result(r)
         import CElementGeometry, ExitStatus
         implicit none
         type(CElementGeometry), intent(in) :: gd
         type(ExitStatus) :: r
       end function print_element_description_wrapper
    end interface
    type(ElementGeometry), intent(in) :: gd
    type(ExitStatus) :: r
    r = print_element_description_wrapper(convert_element_geometry(gd))
  end function print_element_description
  !
  function get_element_description(e, l, f, d) result(r)
    use, intrinsic :: iso_c_binding, only: c_loc
    implicit none
    interface
       function get_element_description_wrapper(e, l, f, d) &
            bind(c,name = 'castem_hho_get_element_description') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_char
         import ElementDescription, CElementGeometry, ExitStatus
         implicit none
         type(ElementDescription), intent(out):: e
         character(len=1,kind=c_char), dimension(*), intent(in) :: l
         character(len=1,kind=c_char), dimension(*), intent(in) :: f
         type(CElementGeometry), intent(in) :: d
         type(ExitStatus) :: r
       end function get_element_description_wrapper
    end interface
    type(ElementDescription), intent(out):: e
    character(len=*), intent(in) :: l
    character(len=*), intent(in) :: f
    type(ElementGeometry), intent(in) :: d
    type(ExitStatus) :: r
    r = get_element_description_wrapper(e, &
         convert_fortran_string(l), &
         convert_fortran_string(f), &
         convert_element_geometry(d))
  end function get_element_description
  !
!  function compute_gradients(g, ed, eg, ewk, swk, ufaces) result(r)
!    use, intrinsic :: iso_c_binding, only: c_loc
!    implicit none
!    interface
!       function compute_gradients_wrapper(g , ed, eg, ewk, swk, ufaces) &
!            bind(c,name = 'castem_hho_compute_gradients') &
!            result(r)
!         use, intrinsic :: iso_c_binding, only: c_ptr
!         import ElementDescription, CElementGeometry, ExitStatus
!         implicit none
!         type(c_ptr), intent(in), value :: g
!         type(ElementDescription), intent(in):: ed
!         type(CElementGeometry), intent(in) :: eg
!         type(c_ptr), intent(in), value :: ewk
!         type(c_ptr), intent(in), value :: swk
!         type(c_ptr), intent(in), value :: ufaces
!         type(ExitStatus) :: r
!       end function compute_gradients_wrapper
!    end interface
!    double precision, dimension (:), target :: g;
!    type(ElementDescription), intent(in) :: ed
!    type(ElementGeometry), intent(in) :: eg
!    double precision, dimension (:), target :: ewk;
!    double precision, dimension (:), target :: swk;
!    double precision, dimension (:), target :: ufaces;
!    type(ExitStatus) :: r
!    r = compute_gradients_wrapper(c_loc(g), ed, &
!         convert_element_geometry(eg), &
!         c_loc(ewk), c_loc(swk), c_loc(ufaces))
!  end function compute_gradients
    !
    function initialize_workspace(ed, ewk, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, ewk, eg) &
                    bind(c,name = 'castem_hho_initialize_workspace') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: ewk
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: ewk;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(ewk), convert_element_geometry(eg))
    end function initialize_workspace
    !
    function get_gauss_data(ed, pts, wts, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, pts, wts, eg) &
                    bind(c,name = 'castem_hho_get_gauss_data') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: pts
                type(c_ptr), intent(in), value :: wts
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: pts;
        double precision, dimension (:), target :: wts;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(pts), c_loc(wts), convert_element_geometry(eg))
    end function get_gauss_data
    !
    function compute_gradients(ed, ewk, g, ufaces, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, ewk, g, ufaces, eg) &
                    bind(c,name = 'castem_hho_compute_gradients') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: ewk
                type(c_ptr), intent(in), value :: g
                type(c_ptr), intent(in), value :: ufaces
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: ewk;
        double precision, dimension (:), target :: g;
        double precision, dimension (:), target :: ufaces;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(ewk), c_loc(g), c_loc(ufaces), convert_element_geometry(eg))
    end function compute_gradients
    !
    function compute_internal_forces(ed, ewk, out, in1, in2, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, ewk, out, in1, in2, eg) &
                    bind(c,name = 'castem_hho_compute_internal_forces') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: ewk
                type(c_ptr), intent(in), value :: out
                type(c_ptr), intent(in), value :: in1
                type(c_ptr), intent(in), value :: in2
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: ewk;
        double precision, dimension (:), target :: out;
        double precision, dimension (:), target :: in1;
        double precision, dimension (:), target :: in2;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(ewk), c_loc(out), c_loc(in1), c_loc(in2), convert_element_geometry(eg))
    end function compute_internal_forces
    !
    function compute_system(ed, ewk, out, in1, in2, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, ewk, out, in1, in2, eg) &
                    bind(c,name = 'castem_hho_compute_system') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: ewk
                type(c_ptr), intent(in), value :: out
                type(c_ptr), intent(in), value :: in1
                type(c_ptr), intent(in), value :: in2
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: ewk;
        double precision, dimension (:), target :: out;
        double precision, dimension (:), target :: in1;
        double precision, dimension (:), target :: in2;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(ewk), c_loc(out), c_loc(in1), c_loc(in2), convert_element_geometry(eg))
    end function compute_system
    !
    function decondensate(ed, ewk, in1, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(ed, ewk, in1, eg) &
                    bind(c,name = 'castem_hho_decondensate') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementDescription, CElementGeometry, ExitStatus
                implicit none
                type(ElementDescription), intent(in):: ed
                type(c_ptr), intent(in), value :: ewk
                type(c_ptr), intent(in), value :: in1
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementDescription), intent(in) :: ed
        double precision, dimension (:), target :: ewk;
        double precision, dimension (:), target :: in1;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(ed, c_loc(ewk), c_loc(in1), convert_element_geometry(eg))
    end function decondensate
  !
end module castem_hho
