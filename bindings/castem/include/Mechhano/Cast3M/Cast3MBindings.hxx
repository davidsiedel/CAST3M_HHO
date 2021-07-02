//
// Created by ds261764 on 27/01/2021.
//

#ifndef MECHHCANO_CASTEM_BINDINGS_HXX
#define MECHHCANO_CASTEM_BINDINGS_HXX

#include "Mechhano/Cast3M/Config.hxx"
#include<cstdint>

extern "C" {

    typedef enum {
      CASTEM_HHO_SUCCESS = 0,
      CASTEM_HHO_FAILURE = 1
    } castem_hho_exit_flag;

    struct castem_hho_status {
      int exit_status;
      char *const msg;
    };

    MECHHANO_CASTEM_EXPORT castem_hho_status castem_hho_report_success();

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_report_failure(const char *const);

    MECHHANO_CASTEM_EXPORT castem_hho_status castem_hho_handle_cxx_exception();

    // forward declaration
    //struct castem_hho_element_description;
    // forward declaration
    //struct castem_hho_element_geometry;

    struct castem_hho_element_geometry {
        //! \brief euclidian dimension
        int64_t dim_eucli;
        //! \brief number of vertices_coordinates of the element
        int64_t num_vertices;
        //! \brief number of vertices_coordinates of the faces
        int64_t num_faces;
        //! \brief number of vertices_coordinates per faces
        int64_t *num_vertices_per_face;
        //! \brief faces connectivity ; f0_node0, ..., fN_node0, ..., f0_node1, ..., fN_node1
        int64_t *connectivity;
        //! \brief coordinates of the vertices_coordinates ; x0, ..., xN, y0, ..., yN
        double  *vertices_coordinates;
    };

    typedef castem_hho_status (*castem_hho_compute_gauss_data_ptr)(
            double *const, // points
            double *const, // weights
            const castem_hho_element_geometry *const
    );

    typedef castem_hho_status (*castem_hho_init_workspace_ptr)(
            double *const, // worksapce
            const castem_hho_element_geometry *const
    );

    typedef castem_hho_status (*castem_hho_compute_gradients_ptr)(
            double *const, // worksapce
            double *const, // OUT : gradients
            const double *const, // IN face_unknown
            const castem_hho_element_geometry *const
    );

    typedef castem_hho_status (*castem_hho_compute_internal_forces_ptr)(
            double *const, // worksapce
            double *const, // OUT
            const double *const, // IN : stress
            const double *const, // IN : face_unknown
            const castem_hho_element_geometry *const
    );

    typedef castem_hho_status (*castem_hho_compute_system_ptr)(
            double *const, // worksapce
            double *const, // OUT
            const double *const, // IN : Ktan
            const double *const, // IN : Res
            const castem_hho_element_geometry *const
    );

    typedef castem_hho_status (*castem_hho_decondensate_ptr)(
            double *const, // worksapce
            const double *const, // IN : faces correction
            const castem_hho_element_geometry *const
    );

    struct castem_hho_element_description {
        int64_t element_workspace_size;
        int64_t num_cell_quadrature_points;
        int64_t element_size;
        int64_t face_size;
        // --- functions
        castem_hho_compute_gauss_data_ptr get_gauss_data;
        castem_hho_init_workspace_ptr initialize_workspace;
        castem_hho_compute_gradients_ptr compute_gradients;
        castem_hho_compute_internal_forces_ptr compute_internal_forces;
        castem_hho_compute_system_ptr compute_system;
        castem_hho_decondensate_ptr decondensate;
    };

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_print_element_geometry(const castem_hho_element_geometry *const);

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_element_description(castem_hho_element_description *const,
                                       const char *const,
                                       const char *const,
                                       const castem_hho_element_geometry *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_initialize_workspace(const castem_hho_element_description *const,
                                    double *const, //worksapce
                                    const castem_hho_element_geometry *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_gauss_data(const castem_hho_element_description *const ed,
                              double *const pts, // points
                                          double *const wts, // weights
                                          const castem_hho_element_geometry *const gd
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_compute_gradients(const castem_hho_element_description *const,
                                 double *const, //worksapce
                                 double *const, // OUT : gradients
                                 const double *const, // IN : face_unknown
                                 const castem_hho_element_geometry *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_compute_internal_forces(const castem_hho_element_description *const,
                                       double *const, //worksapce
                                       double *const, // OUT
                                       const double *const, // IN : stress
                                       const double *const, // IN : face_unknown
                                       const castem_hho_element_geometry *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_compute_system(const castem_hho_element_description *const,
                              double *const, //worksapce
                              double *const, // OUT
                              const double *const, // IN : Ktan
                              const double *const, // IN : Res
                              const castem_hho_element_geometry *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_decondensate(const castem_hho_element_description *const,
                              double *const, //worksapce
                              const double *const, // IN : faces correction
                              const castem_hho_element_geometry *const
    );

}  // end of extern "C"

//    /*!
//     * \brief compute gradient
//     * \param[in] pointer to a
//     * \param[in] pointer to a
//     */
//    typedef castem_hho_status (*castem_hho_compute_gradients_ptr)(
//            double *const, //worksapce
//            double *const, // OUT : gradient data
//            const double *const, // IN : face unknowns data
//            const castem_hho_element_geometry *const
//    );
//
//    /*!
//     * \brief compute internal forces
//     * \param[in] pointer to a
//     * \param[in] pointer to a
//     */
//    typedef castem_hho_status (*castem_hho_compute_internal_forces_ptr)(
//            double *const, //worksapce
//            double *const, // OUT : stress data
//            const double *const, // IN : face unknowns data
//            const castem_hho_element_geometry *const
//    );
//
//    /*!
//     * \brief compute element system
//     * \param[in] pointer to a
//     * \param[in] pointer to a
//     */
//    typedef castem_hho_status (*castem_hho_compute_system_ptr)(
//            double *const, //worksapce
//            double *const, // OUT : system data
//            const double *const, // IN : Ktan data, internal and external forces data
//            const castem_hho_element_geometry *const
//    );

//    struct castem_hho_computation_data{
//        //! \brief stabilization parameter
//        double stabilization_parameter;
//        //! \brief face displacement values
//        double* face_displacement_values;
//        //! \brief grad values ; g00, g01, g02, g10, g11, g12, g21, g20, g21, g22
//        double* gradient_values;
//        //! \brief stress values ; s00, s11, s22, s01, s02, s12
//        double* stress_values;
//        //! \brief Ktan values ; row-wise, 6 * 6 values
//        double* tangent_operator_values;
//        //! \brief element tangent matrix ; row-wise, element_size * element_size values
//        double* element_matrix_values;
//        //! \brief element tangent residual ; row-wise, element_size values
//        double* element_residual_values;
//    };

//    struct castem_hho_element_workspace{
//        int64_t precompute;
//        double* workspace;
//    };

//    struct castem_hho_integration_point_data {
//        //! \brief number of cell quadrature points
//        int64_t num_cell_quadrature_points;
//        //! \brief cell quadrature points coordinates ; x0, ..., xN, y0, ..., yN
//        double* cell_quadrature_points;
//        //! \brief cell quadrature weights ; w0, ..., wN
//        double* cell_quadrature_weights;
//    };

    /*!
     * \brief fill the worskapce array with relevant values (B matrices values, cell unknowns, etc.)
     * \param[in] pointer to a
     * \param[in] pointer to a
     */
//    typedef castem_hho_status (*castem_hho_initialize_workspace_ptr)(
//            double *const,
//            const castem_hho_element_geometry *const
//            );



#endif  // MECHHCANO_CASTEM_BINDINGS_HXX
