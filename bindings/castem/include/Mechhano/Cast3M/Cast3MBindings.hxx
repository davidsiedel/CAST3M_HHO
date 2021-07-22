//
// Created by ds261764 on 27/01/2021.
//

#ifndef MECHHCANO_CASTEM_BINDINGS_HXX
#define MECHHCANO_CASTEM_BINDINGS_HXX
#define nfhmax 20

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
        //! \brief faces connectivity ; f0_node0, ..., fN_node0, ..., f0_node1, ..., fN_node1
//        int64_t *connectivity;
        //! \brief coordinates of the vertices_coordinates ; x0, ..., xN, y0, ..., yN
//        double  *vertices_coordinates;
        int64_t connectivity [nfhmax * 2];
        double vertices_coordinates [nfhmax * 2];
    };


    /*!
    * \brief generic element description
    *
    */
    struct castem_hho_element_description {
        int64_t dim_eucli;
        int64_t num_vertices;
        int64_t num_faces;
        int64_t num_vertices_per_face [nfhmax];
        int64_t num_quadrature_points;
        int64_t dir_dof_face_unknown;
        int64_t dir_dof_cell_unknown;
        int64_t dir_dof_gradient;
        int64_t dim_field;
        int64_t dim_space_element;
        int64_t dim_face_block;
        int64_t dim_cell_block;
        int64_t dim_MB;
        int64_t dim_MB_matrices;
        int64_t dim_MSTAB;
        int64_t dim_MKCC;
        int64_t dim_MKCF;
        int64_t dim_MVC;
    };

    typedef castem_hho_status (*castem_hho_get_gradient_operator_ptr)(
            const castem_hho_element_geometry *const,
            double *const, // data
            int64_t // index
    );

    typedef castem_hho_status (*castem_hho_get_stabilization_operator_ptr)(
            const castem_hho_element_geometry *const,
            double *const // data
    );

    typedef castem_hho_status (*castem_hho_get_gauss_weight_ptr)(
            const castem_hho_element_geometry *const,
            double *const, // data
            int64_t // index
    );

    typedef castem_hho_status (*castem_hho_get_gauss_point_ptr)(
            const castem_hho_element_geometry *const,
            double *const, // data
            int64_t // index
    );

    typedef castem_hho_status (*castem_hho_invert_matrix_ptr)(
            double *const, // data
            int64_t // index
    );

    struct castem_hho_element_functions {
        castem_hho_get_gradient_operator_ptr get_gradient_operator;
        castem_hho_get_stabilization_operator_ptr get_stabilization_operator;
        castem_hho_get_gauss_weight_ptr get_gauss_weight;
        castem_hho_get_gauss_point_ptr get_gauss_point;
    };

    struct castem_hho_generic_functions {
        castem_hho_invert_matrix_ptr invert_matrix;
    };

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_element_description(castem_hho_element_description *const,
                                       const char *const,
                                       const char *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_element_functions(castem_hho_element_functions *const,
//                                     const castem_hho_element_geometry *const,
                                       const char *const,
                                       const char *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_generic_functions(castem_hho_generic_functions *const,
    //                                     const castem_hho_element_geometry *const,
                                     const char *const,
                                     const char *const
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_gradient_operator(const castem_hho_element_functions *const,
                                     const castem_hho_element_geometry *const,
                                     double *const, // data
                                     int64_t // index
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_stabilization_operator(const castem_hho_element_functions *const,
                                     const castem_hho_element_geometry *const,
                                     double *const // data
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_gauss_weight(const castem_hho_element_functions *const,
                                     const castem_hho_element_geometry *const,
                                    double *const, // data
                                     int64_t // index
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_get_gauss_point(const castem_hho_element_functions *const,
                                     const castem_hho_element_geometry *const,
                                     double *const, // data
                                     int64_t // index
    );

    MECHHANO_CASTEM_EXPORT castem_hho_status
    castem_hho_invert_matrix(
//            const castem_hho_element_functions *const,
                                const castem_hho_generic_functions *const,
                                double *const, // data
                                int64_t // index
    );

}  // end of extern "C"

#endif  // MECHHCANO_CASTEM_BINDINGS_HXX
