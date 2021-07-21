module castem_hho
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int, c_ptr, c_funptr, c_char, c_null_ptr, c_loc
    use, intrinsic :: iso_fortran_env, only: int64
    use castem_hho_fortran_utilities
    !
    !! Maximum Number of Vertices for a face (3D)/an edge (2D) of a HHO cell
    !! (must be consistent with the same variable in C3m include CCHHOPR.INC)
    integer(c_int64_t), parameter :: nfhmax = 20
    !
    type, bind(c) :: ElementGeometry
        integer(c_int64_t), dimension(nfhmax * 2) :: connectivity
        double precision, dimension(nfhmax * 2) :: vertices_coordinates
    end type ElementGeometry
    !
    type, bind(c) :: ElementDescription
        integer(c_int64_t) :: dim_eucli
        integer(c_int64_t) :: num_vertices
        integer(c_int64_t) :: num_faces
        integer(c_int64_t), dimension(nfhmax) :: num_vertices_per_face
        integer(c_int64_t) :: num_quadrature_points
        integer(c_int64_t) :: dir_dof_face_unknown
        integer(c_int64_t) :: dir_dof_cell_unknown
        integer(c_int64_t) :: dir_dof_gradient
        !!        integer(c_int64_t) :: dim_gradient     !! Always 9
        integer(c_int64_t) :: dim_field          !! = dim_eucli for mechanics (= 1 for phase-field or thermics)
        integer(c_int64_t) :: dim_space_element  !! dim_field * (dir_dof_face_unknown * num_faces + dir_dof_cell_unknown)
        integer(c_int64_t) :: dim_face_block     !! dim_field * (dir_dof_face_unknown * num_faces)
        integer(c_int64_t) :: dim_cell_block     !! dim_field * (dir_dof_cell_unknown)
        integer(c_int64_t) :: dim_MB             !! dim.Matri9, dim_space_element)
        integer(c_int64_t) :: dim_MB_matrices    !! dim_MB * num_quadrature_points
        integer(c_int64_t) :: dim_MSTAB          !! dim.matridim_space_element,dim_space_element)
        integer(c_int64_t) :: dim_MKCC           !! dim.Matridim_cell_block,dim_cell_block)
        integer(c_int64_t) :: dim_MKCF           !! dim.Matridim_cell_block,dim_face_block)
        integer(c_int64_t) :: dim_MVC            !! dim.Vector(dim_cell_block) [Vector(.) = Matri1,.)]
    end type ElementDescription
    !
    type, bind(c) :: ElementFunctions
        type(c_funptr) get_gradient_operator
        type(c_funptr) get_stabilization_operator
        type(c_funptr) get_gauss_weight
        type(c_funptr) get_gauss_point
    end type ElementFunctions
    !
    type, bind(c) :: GenericFunctions
        type(c_funptr) invert_matrix
    end type GenericFunctions
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
    ! ------------------------------------------------------------------------------------------------------------------
    ! BUILD ELEEMNT DESCRIPTION OBJECT
    ! ------------------------------------------------------------------------------------------------------------------
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
    ! ------------------------------------------------------------------------------------------------------------------
    ! BUILD ELEMENT FUNCTIONS OBJECT
    ! ------------------------------------------------------------------------------------------------------------------
    function get_element_functions(elem_funs, l, f) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function get_element_functions_wrapper(elem_funs, l, f) &
                    bind(c,name = 'castem_hho_get_element_functions') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_char
                import ElementFunctions, ExitStatus
                implicit none
                type(ElementFunctions), intent(out):: elem_funs
                character(len=1,kind=c_char), dimension(*), intent(in) :: l
                character(len=1,kind=c_char), dimension(*), intent(in) :: f
                type(ExitStatus) :: r
            end function get_element_functions_wrapper
        end interface
        type(ElementFunctions), intent(out):: elem_funs
        character(len=*), intent(in) :: l
        character(len=*), intent(in) :: f
        type(ExitStatus) :: r
        r = get_element_functions_wrapper(elem_funs, &
                convert_fortran_string(l), &
                convert_fortran_string(f))
    end function get_element_functions
    ! ------------------------------------------------------------------------------------------------------------------
    ! BUILD GENERIC FUNCTIONS OBJECT
    ! ------------------------------------------------------------------------------------------------------------------
    function get_generic_functions(gene_funs, l, f) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function get_generic_functions_wrapper(gene_funs, l, f) &
                    bind(c,name = 'castem_hho_get_generic_functions') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_char
                import GenericFunctions, ExitStatus
                implicit none
                type(GenericFunctions), intent(out):: gene_funs
                character(len=1,kind=c_char), dimension(*), intent(in) :: l
                character(len=1,kind=c_char), dimension(*), intent(in) :: f
                type(ExitStatus) :: r
            end function get_generic_functions_wrapper
        end interface
        type(GenericFunctions), intent(out):: gene_funs
        character(len=*), intent(in) :: l
        character(len=*), intent(in) :: f
        type(ExitStatus) :: r
        r = get_generic_functions_wrapper(gene_funs, &
                convert_fortran_string(l), &
                convert_fortran_string(f))
    end function get_generic_functions
    ! ------------------------------------------------------------------------------------------------------------------
    ! COMPUTE GRADIENT OPERATOR AT A SINGLE GAUSS POINT (INDEXED I)
    ! ------------------------------------------------------------------------------------------------------------------
    function get_gradient_operator(elem_funs, elem_geom, data, index) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc, c_int64_t
        implicit none
        interface
            function get_gradient_operator_wrapper(elem_funs, elem_geom, data, index) &
                    bind(c,name = 'castem_hho_get_gradient_operator') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t
                import ElementFunctions, ExitStatus, ElementGeometry
                implicit none
                type(ElementFunctions), intent(in) :: elem_funs
                type(ElementGeometry), intent(in) :: elem_geom
                type(c_ptr), intent(in), value :: data
                integer(c_int64_t), intent(in), value :: index
                type(ExitStatus) :: r
            end function get_gradient_operator_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: data
        type(ElementGeometry), intent(in) :: elem_geom
        integer(c_int64_t), intent(in) :: index
        type(ExitStatus) :: r
        r = get_gradient_operator_wrapper(elem_funs, elem_geom, c_loc(data), index)
    end function get_gradient_operator
    ! ------------------------------------------------------------------------------------------------------------------
    ! COMPUTE GAUSS WEIGHT AT A SINGLE GAUSS POINT (INDEXED I)
    ! ------------------------------------------------------------------------------------------------------------------
    function get_gauss_weight(elem_funs, elem_geom, data, index) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc, c_int64_t
        implicit none
        interface
            function get_gauss_weight_wrapper(elem_funs, elem_geom, data, index) &
                    bind(c,name = 'castem_hho_get_gauss_weight') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t
                import ElementFunctions, ExitStatus, ElementGeometry
                implicit none
                type(ElementFunctions), intent(in) :: elem_funs
                type(ElementGeometry), intent(in) :: elem_geom
                type(c_ptr), intent(in), value :: data
                integer(c_int64_t), intent(in), value :: index
                type(ExitStatus) :: r
            end function get_gauss_weight_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: data;
        type(ElementGeometry), intent(in) :: elem_geom
        integer(c_int64_t), intent(in) :: index
        type(ExitStatus) :: r
        r = get_gauss_weight_wrapper(elem_funs, elem_geom, c_loc(data), index)
    end function get_gauss_weight
    ! ------------------------------------------------------------------------------------------------------------------
    ! COMPUTE GAUSS POINT AT A SINGLE GAUSS POINT (INDEXED I)
    ! ------------------------------------------------------------------------------------------------------------------
    function get_gauss_point(elem_funs, elem_geom, data, index) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc, c_int64_t
        implicit none
        interface
            function get_gauss_point_wrapper(elem_funs, elem_geom, data, index) &
                    bind(c,name = 'castem_hho_get_gauss_point') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t
                import ElementFunctions, ExitStatus, ElementGeometry
                implicit none
                type(ElementFunctions), intent(in) :: elem_funs
                type(ElementGeometry), intent(in) :: elem_geom
                type(c_ptr), intent(in), value :: data
                integer(c_int64_t), intent(in), value :: index
                type(ExitStatus) :: r
            end function get_gauss_point_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: data;
        type(ElementGeometry), intent(in) :: elem_geom
        integer(c_int64_t), intent(in) :: index
        type(ExitStatus) :: r
        r = get_gauss_point_wrapper(elem_funs, elem_geom, c_loc(data), index)
    end function get_gauss_point
    ! ------------------------------------------------------------------------------------------------------------------
    ! COMPUTE STABILIZATION OPERATOR
    ! ------------------------------------------------------------------------------------------------------------------
    function get_stabilization_operator(elem_funs, elem_geom, data) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc
        implicit none
        interface
            function get_stabilization_operator_wrapper(elem_funs, elem_geom, data) &
                    bind(c,name = 'castem_hho_get_stabilization_operator') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr
                import ElementFunctions, ExitStatus, ElementGeometry
                implicit none
                type(ElementFunctions), intent(in) :: elem_funs
                type(ElementGeometry), intent(in) :: elem_geom
                type(c_ptr), intent(in), value :: data
                type(ExitStatus) :: r
            end function get_stabilization_operator_wrapper
        end interface
        type(ElementFunctions), intent(in) :: elem_funs
        double precision, dimension (:), target :: data
        type(ElementGeometry), intent(in) :: elem_geom
        type(ExitStatus) :: r
        r = get_stabilization_operator_wrapper(elem_funs, elem_geom, c_loc(data))
    end function get_stabilization_operator
    ! ------------------------------------------------------------------------------------------------------------------
    ! INVERT MATRIX
    ! ------------------------------------------------------------------------------------------------------------------
    function invert_matrix(gene_funs, data, index) result(r)
        use, intrinsic :: iso_c_binding, only: c_loc, c_int64_t
        implicit none
        interface
            function invert_matrix_wrapper(gene_funs, data, index) &
                    bind(c,name = 'castem_hho_invert_matrix') &
                            result(r)
                use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t
                import ExitStatus, GenericFunctions
                implicit none
                type(GenericFunctions), intent(in) :: gene_funs
                type(c_ptr), intent(in), value :: data
                integer(c_int64_t), intent(in), value :: index
                type(ExitStatus) :: r
            end function invert_matrix_wrapper
        end interface
        type(GenericFunctions), intent(in) :: gene_funs
        double precision, dimension (:), target :: data;
        integer(c_int64_t), intent(in) :: index
        type(ExitStatus) :: r
        r = invert_matrix_wrapper(gene_funs, c_loc(data), index)
    end function invert_matrix
end module castem_hho