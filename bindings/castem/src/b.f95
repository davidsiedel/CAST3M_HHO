module castem_hho
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_ptr, c_funptr, c_char, c_null_ptr, c_loc
    use, intrinsic :: iso_fortran_env, only: int64
    use castem_hho_fortran_utilities
    !
    type :: ElementGeometry
        integer(c_int64_t), dimension (:), pointer :: connectivity
        double precision, dimension (:), pointer :: vertices_coordinates
!        integer(c_int64_t), dimension (40) :: connectivity
!        double precision, dimension (40) :: vertices_coordinates
    end type ElementGeometry
    !
    type, bind(c) :: ElementDescription
        integer(c_int64_t) :: dim_eucli
        integer(c_int64_t) :: num_vertices
        integer(c_int64_t) :: num_faces
        integer(c_int64_t), dimension (20) :: num_vertices_per_face
        integer(c_int64_t) :: num_quadrature_points
        integer(c_int64_t) :: dim_space_face_unknown
        integer(c_int64_t) :: dim_space_cell_unknown
        integer(c_int64_t) :: dim_space_gradient
        integer(c_int64_t) :: dim_field
        integer(c_int64_t) :: dim_space_element
        integer(c_int64_t) :: dim_face_block
        integer(c_int64_t) :: dim_cell_block
        integer(c_int64_t) :: dim_B
        integer(c_int64_t) :: dim_B_matrices
        integer(c_int64_t) :: dim_Z
        integer(c_int64_t) :: dim_MCC
        integer(c_int64_t) :: dim_MCF
        integer(c_int64_t) :: dim_VC
    end type ElementDescription
    !
    type, bind(c) :: ElementFunctions
        type(c_funptr) get_gauss_data
        type(c_funptr) compute_operators
        type(c_funptr) compute_gradients
        type(c_funptr) compute_internal_forces
        type(c_funptr) condensate_system
        type(c_funptr) decondensate_system
    end type ElementFunctions
    !
    type, bind(c) :: CElementGeometry
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
    function get_element_description(elem_desc, l, f) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function get_element_description_wrapper(elem_desc, l, f) &
                    bind(c,name = 'castem_hho_get_element_description') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_char
                import ElementDescription, ExitStatus
                implicit none
                type(ElementDescription), intent(out):: elem_desc
                character(len=1,kind=c_char), dimension(*), intent(in) :: l
                character(len=1,kind=c_char), dimension(*), intent(in) :: f
                type(ExitStatus) :: r
            end function get_element_description_wrapper
        end interface
        type(ElementDescription), intent(out):: elem_desc
        character(len=*), intent(in) :: l
        character(len=*), intent(in) :: f
        type(ExitStatus) :: r
        r = get_element_description_wrapper(elem_desc, &
                convert_fortran_string(l), &
                convert_fortran_string(f))
    end function get_element_description
    !
    function get_element_functions(elem_funs, l, f, d) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function get_element_functions_wrapper(elem_funs, l, f, d) &
                    bind(c,name = 'castem_hho_get_element_functions') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_char
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(out):: elem_funs
                character(len=1,kind=c_char), dimension(*), intent(in) :: l
                character(len=1,kind=c_char), dimension(*), intent(in) :: f
                type(CElementGeometry), intent(in) :: d
                type(ExitStatus) :: r
            end function get_element_functions_wrapper
        end interface
        type(ElementFunctions), intent(out):: elem_funs
        character(len=*), intent(in) :: l
        character(len=*), intent(in) :: f
        type(ElementGeometry), intent(in) :: d
        type(ExitStatus) :: r
        r = get_element_functions_wrapper(elem_funs, &
                convert_fortran_string(l), &
                convert_fortran_string(f), &
                convert_element_geometry(d))
    end function get_element_functions
    !
    function compute_operators(elem_funs, ewk, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, ewk, eg) &
                    bind(c,name = 'castem_hho_compute_operators') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: ewk
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: ewk;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(ewk), convert_element_geometry(eg))
    end function compute_operators
    !
    function get_gauss_data(elem_funs, pts, wts, eg) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, pts, wts, eg) &
                    bind(c,name = 'castem_hho_get_gauss_data') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: pts
                type(c_ptr), intent(in), value :: wts
                type(CElementGeometry), intent(in) :: eg
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: pts;
        double precision, dimension (:), target :: wts;
        type(ElementGeometry), intent(in) :: eg
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(pts), c_loc(wts), convert_element_geometry(eg))
    end function get_gauss_data
    !
    function compute_gradients(elem_funs, ops, g, ufaces) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, ops, g, ufaces) &
                    bind(c,name = 'castem_hho_compute_gradients') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: ops
                type(c_ptr), intent(in), value :: g
                type(c_ptr), intent(in), value :: ufaces
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: ops;
        double precision, dimension (:), target :: g;
        double precision, dimension (:), target :: ufaces;
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(ops), c_loc(g), c_loc(ufaces))
    end function compute_gradients
    !
    function compute_internal_forces(elem_funs, ops, out, in1, in2) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, ops, out, in1, in2) &
                    bind(c,name = 'castem_hho_compute_internal_forces') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: ops
                type(c_ptr), intent(in), value :: out
                type(c_ptr), intent(in), value :: in1
                type(c_ptr), intent(in), value :: in2
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: ops;
        double precision, dimension (:), target :: out;
        double precision, dimension (:), target :: in1;
        double precision, dimension (:), target :: in2;
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(ops), c_loc(out), c_loc(in1), c_loc(in2))
    end function compute_internal_forces
    !
    function compute_system(elem_funs, ops, out, in1, in2) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, ops, out, in1, in2) &
                    bind(c,name = 'castem_hho_compute_system') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: ops
                type(c_ptr), intent(in), value :: out
                type(c_ptr), intent(in), value :: in1
                type(c_ptr), intent(in), value :: in2
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: ops;
        double precision, dimension (:), target :: out;
        double precision, dimension (:), target :: in1;
        double precision, dimension (:), target :: in2;
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(ops), c_loc(out), c_loc(in1), c_loc(in2))
    end function compute_system
    !
    function decondensate(elem_funs, ops, in1) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function compute_gradients_wrapper(elem_funs, ops, in1) &
                    bind(c,name = 'castem_hho_decondensate') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, CElementGeometry, ExitStatus
                implicit none
                type(ElementFunctions), intent(in):: elem_funs
                type(c_ptr), intent(in), value :: ops
                type(c_ptr), intent(in), value :: in1
                type(ExitStatus) :: r
            end function compute_gradients_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: ops;
        double precision, dimension (:), target :: in1;
        type(ExitStatus) :: r
        r = compute_gradients_wrapper(elem_funs, c_loc(ops), c_loc(in1))
    end function decondensate
    !
end module castem_hho
