//
// Created by dsiedel on 01/07/2021.
//

#ifndef FETA_FUNCTIONS_HXX
#define FETA_FUNCTIONS_HXX

#define verbose false

#include "Mechhano/Cast3M/Cast3MBindings.hxx"
#include "feta/hybrid.hxx"

namespace feta::hybrid{

    castem_hho_status invertMatrix(
            double *const data, // worksapce
            int64_t index
    ) {
        try {
            if (verbose) std::cout << "******** INVERTING MATRIX" << std::endl;
            Eigen::Matrix<real , Eigen::Dynamic, Eigen::Dynamic> mat(index, index);
            for (int i = 0; i < index; ++i) {
                for (int j = 0; j < index; ++j) {
                    mat(i,j) = data[j + index * i];
                }
            }
            Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> eye = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Identity(index, index);
            Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> inv = mat.llt().solve(eye);
            for (int i = 0; i < index; ++i) {
                for (int j = 0; j < index; ++j) {
                    data[j + index * i] = inv(i, j);
                }
            }
            if (verbose) std::cout << "******** MATRIX INVERTED" << std::endl;
        } catch (...) {
            return castem_hho_handle_cxx_exception();
        }
        return castem_hho_report_success();
    }

    static castem_hho_status getGenericFunctions(
            castem_hho_generic_functions *const gene_funs
    ){
        try {
            if (verbose) std::cout << "******** INITIALIZING GENERIC FUNCTIONS" << std::endl;
            gene_funs->invert_matrix = invertMatrix;
            if (verbose) std::cout << "******** GENERIC FUNCTIONS INITIALIZED" << std::endl;
        } catch (...) {
            return castem_hho_handle_cxx_exception();
        }
        return castem_hho_report_success();
    }

//    template<ElementType elem_t, FieldType field_t, intg k, intg l>
    template<ElementType elem_t, FieldType field_t, intg l, intg k>
    struct Cas3MFunctions{
        using Cast3MElement = Cast3MElem<elem_t, field_t, k, l>;
        static constexpr intg num_vertices = Elem2<elem_t, k>::num_nodes;
        static constexpr intg grad_size = Cast3MElement::grad_size;
        static constexpr intg elem_size = Cast3MElement::elem_size;
        static constexpr intg nq = Cast3MElement::CellCmpt::dim_q;
        static constexpr intg cell_block_size = Cast3MElement::Space ::cl * Cast3MElement::field_size;
        static constexpr intg face_block_size = Cast3MElement::elem_size - cell_block_size;

        static std::array<intg, nfhmax> getNumVerticesPerFace() {
            std::array<intg, nfhmax> num_vertices_per_face;
            for (int i = 0; i < nfhmax; ++i) {
                if (i < num_vertices) {
                    num_vertices_per_face[i] = 2;
                }
                else num_vertices_per_face[i] = 0;
            }
            return num_vertices_per_face;
        }

        static castem_hho_status getDescription(
                castem_hho_element_description *const d
        ){
            try {
                if (verbose) std::cout << "******** INITIALIZING ELEMENT DESCRIPTION" << std::endl;
//                d->dim_eucli = Cas3MFunctions::dim_eucli;
                d->dim_eucli = Cast3MElement::d;
                d->num_vertices = Elem2<elem_t, k>::num_nodes;
                d->num_faces = Elem2<elem_t, k>::num_nodes;
                std::array<intg, nfhmax> tab = getNumVerticesPerFace();
                for (int i = 0; i < nfhmax; ++i) {
                    d->num_vertices_per_face[i] = tab[i];
                }
                d->num_quadrature_points = Cast3MElement::CellCmpt::dim_q;
                d->dir_dof_face_unknown = Cast3MElement::fk;
                d->dir_dof_cell_unknown = Cast3MElement::cl;
                d->dir_dof_gradient = Cast3MElement::ck;
                d->dim_field = Cast3MElement::field_size;
                d->dim_space_element = Cast3MElement::elem_size;
                d->dim_face_block = d->dim_field * d->num_faces * d->dir_dof_face_unknown;
                d->dim_cell_block = d->dim_field * d->dir_dof_cell_unknown;
                d->dim_MB = d->dim_space_element * Cast3MElement::grad_size;
                d->dim_MB_matrices = d->dim_space_element * Cast3MElement::grad_size * d->num_quadrature_points;
                d->dim_MSTAB = d->dim_space_element * d->dim_space_element;
                d->dim_MKCC = d->dim_cell_block * d->dim_cell_block;
                d->dim_MKCF = d->dim_cell_block * d->dim_face_block;
                d->dim_MVC = d->dim_cell_block;
                if (verbose) std::cout << "******** ELEMENT DESCRIPTION INITIALIZED" << std::endl;
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_gradient_operator(
                const castem_hho_element_geometry *const gd,
                double *const data, // worksapce
                int64_t index
        ) {
            try {
                if (verbose) std::cout << "******** GETTING GRADIENT OPERATOR AT GAUSS POINT : " << index << std::endl;
                // --- BUILDING THE ELEMENT
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
                if (verbose) std::cout << "- VERICES AS EXTRACTED :" << std::endl;
                if (verbose) std::cout << verts << std::endl;
                if (verbose) std::cout << "- CONNECTIVITY AS EXTRACTED :" << std::endl;
                if (verbose) std::cout << ordening << std::endl;
                Cast3MElement elem = Cast3MElement(verts, ordening);
                // --- FILLING THE ELEMENT WORKSPACE FOR GRADIENTS
                int effetcive_index = index - 1;
                for (int i = 0; i < grad_size; ++i) {
                    for (int j = 0; j < elem_size; ++j) {
                        int array_index = i * elem_size + j;
                        data[array_index] = elem.gradient_operators[effetcive_index](i, j);
                    }
                }
                if (verbose) std::cout << elem.gradient_operators[effetcive_index] << std::endl;
                if (verbose) std::cout << "******** GRADIENT OPERATOR AT GAUSS POINT : " << index << " COMPUTED" << std::endl;
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_stabilization_operator(
                const castem_hho_element_geometry *const gd,
                double *const data // worksapce
        ) {
            try {
                if (verbose) std::cout << "******** GETTING STABILIZATION OPERATOR" << std::endl;
                // --- BUILDING THE ELEMENT
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
                if (verbose) std::cout << "- VERICES AS EXTRACTED :" << std::endl;
                if (verbose) std::cout << verts << std::endl;
                if (verbose) std::cout << "- CONNECTIVITY AS EXTRACTED :" << std::endl;
                if (verbose) std::cout << ordening << std::endl;
                Cast3MElement elem = Cast3MElement(verts, ordening);
                // --- FILLING THE ELEMENT WORKSPACE FOR GRADIENTS
                for (int i = 0; i < elem_size; ++i) {
                    for (int j = 0; j < elem_size; ++j) {
                        int array_index = i * elem_size + j;
                        data[array_index] = elem.stabilization_operator(i, j);
                    }
                }
                if (verbose) std::cout << elem.stabilization_operator << std::endl;
                if (verbose) std::cout << "******** STABILIZATION OPERATOR COMPUTED" << std::endl;
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_gauss_weight(
                const castem_hho_element_geometry *const gd,
                double *const data, // worksapce
                int64_t index
        ) {
            try {
                if (verbose) std::cout << "******** GETTING GAUSS WEIGHT : " << index << std::endl;
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
                Cast3MElement elem = Cast3MElement(verts, ordening);
                int effetcive_index = index - 1;
                real qw = elem.cell_cmpt.getQuadWeight(effetcive_index);
                data[0] = qw;
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_gauss_point(
                const castem_hho_element_geometry *const gd,
                double *const data, // worksapce
                int64_t index
        ) {
            try {
                if (verbose) std::cout << "******** GETTING GAUSS POINT : " << index << std::endl;
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
                Cast3MElement elem = Cast3MElement(verts, ordening);
                int effetcive_index = index - 1;
                EigMat<2, 1> qp = elem.cell_cmpt.getQuadPoint(effetcive_index);
                for (int i = 0; i < Cast3MElement::d; ++i) {
                    data[i] = qp(i);
                }
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_nodal_cell_displacement(
                const castem_hho_element_geometry *const gd,
                double *const data, // data
                double *const unknowns // data
        ) {
            try {
                if (verbose) std::cout << "******** GETTING NODAL DISPLACEMENT : " << std::endl;
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
                const intg cell_field_dim = Cast3MElement::cl * Cast3MElement::field_size;
                const intg nodal_displacement_dim = Cast3MElement::num_nodes * Cast3MElement::field_size;
                if (verbose) std::cout << "******** cell_field_dim : " << cell_field_dim << std::endl;
                const EigMapC<cell_field_dim, 1> Uvect(unknowns);
                Cast3MElement elem = Cast3MElement(verts, ordening);
                auto c = elem.cell.getBarycenter();
                auto b = elem.cell.getBounds();
                using Interpolator = EigMat<nodal_displacement_dim, cell_field_dim>;
                Interpolator displacement_interpolation_operator = Interpolator::Zero();
                for (int i = 0; i < Cast3MElement::num_nodes; ++i) {
                    EigMat<Cast3MElement::cl, 1> vectt = elem.space.cell_basis.GetEvaluationVector(verts.col(i), c, b);
                    for (int j = 0; j < Cast3MElement::field_size; ++j) {
                        int rowpos = i + j * Cast3MElement::field_size;
                        int colpos = j * Cast3MElement::cl;
                        displacement_interpolation_operator.template block<1, Cast3MElement::cl>(rowpos, colpos) = vectt;
                    }
                }
                EigMat<Cast3MElement::num_nodes * Cast3MElement::field_size, 1> Dvect = displacement_interpolation_operator * Uvect;
                for (int i = 0; i < Cast3MElement::num_nodes * Cast3MElement::field_size; ++i) {
                    data[i] = Dvect(i);
                }
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status get_nodal_faces_displacement(
                const castem_hho_element_geometry *const gd,
                double *const data, // data
                double *const unknowns // data
        ) {
            try {
                if (verbose) std::cout << "******** GETTING NODAL DISPLACEMENT : " << std::endl;
                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
                const EigMapIntC<2, Cast3MElement::num_nodes> ordening = EigMapIntC<2, Cast3MElement::num_nodes>(gd->connectivity);
//                const EigMap<2, Cast3MElement::num_nodes> ordening = EigMap<2, Cast3MElement::num_nodes>(gd->connectivity);
                const intg cell_field_dim = Cast3MElement::cl * Cast3MElement::field_size;
                const intg nodal_displacement_dim = Cast3MElement::num_nodes * Cast3MElement::field_size;
                if (verbose) std::cout << "******** cell_field_dim : " << cell_field_dim << std::endl;
                const EigMapC<cell_field_dim, 1> Uvect(unknowns);
                Cast3MElement elem = Cast3MElement(verts, ordening);
                auto c = elem.cell.getBarycenter();
                auto b = elem.cell.getBounds();
                using Interpolator = EigMat<nodal_displacement_dim, cell_field_dim>;
                Interpolator displacement_interpolation_operator = Interpolator::Zero();
                for (int i = 0; i < Cast3MElement::num_nodes; ++i) {
                    EigMat<Cast3MElement::cl, 1> vectt = elem.space.cell_basis.GetEvaluationVector(verts.col(i), c, b);
                    for (int j = 0; j < Cast3MElement::field_size; ++j) {
                        int rowpos = i + j * Cast3MElement::field_size;
                        int colpos = j * Cast3MElement::cl;
                        displacement_interpolation_operator.template block<1, Cast3MElement::cl>(rowpos, colpos) = vectt;
                    }
                }
                EigMat<Cast3MElement::num_nodes * Cast3MElement::field_size, 1> Dvect = displacement_interpolation_operator * Uvect;
                for (int i = 0; i < Cast3MElement::num_nodes * Cast3MElement::field_size; ++i) {
                    data[i] = Dvect(i);
                }
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

        static castem_hho_status getFunctions(
                castem_hho_element_functions *const elem_funs
        ){
            try {
                if (verbose) std::cout << "******** INITIALIZING ELEMENT FUNCTIONS" << std::endl;
                elem_funs->get_gradient_operator = get_gradient_operator;
                elem_funs->get_stabilization_operator = get_stabilization_operator;
                elem_funs->get_gauss_weight = get_gauss_weight;
                elem_funs->get_gauss_point = get_gauss_point;
                elem_funs->get_nodal_cell_displacement = get_nodal_cell_displacement;
                if (verbose) std::cout << "******** ELEMENT FUNCTIONS INITIALIZED" << std::endl;
            } catch (...) {
                return castem_hho_handle_cxx_exception();
            }
            return castem_hho_report_success();
        }

//    };

//        static intg getElementOperatorsSize() {
//            // --- GET B MATRICES SIZE
//            const intg B_matrices_size = grad_size * elem_size * nq;
//            // --- GET STABILIZATION MATRIX SIZE
//            const intg Z_matrix_size = elem_size * elem_size;
//            // --- GET TOTAL OPERATORS SIZE
//            const intg elem_ops_size = B_matrices_size + Z_matrix_size;
//            return elem_ops_size;
//        }



//        static intg getWorkSpaceSize() {
//            // --- DEFINE TEMPLATED STRUCTS
//            // --- DETERMINE THE WHOLE ELEMENT WORKSPACE SIZE
//            int element_worskpace_size = 1;
//            // --- operators size
//            element_worskpace_size += getElementOperatorsSize();
//            // --- cell-cell matrix invert size
//            element_worskpace_size += cell_block_size * cell_block_size;
//            // --- cell-faces matrix size
//            element_worskpace_size += cell_block_size * face_block_size;
//            // --- cell vector size
//            element_worskpace_size += cell_block_size;
//            // --- cell unknowns
//            element_worskpace_size += cell_block_size;
//            // --- face lagrange flags
////  element_worskpace_size += n_fcs * Field::field_size;
//            return element_worskpace_size;
//        }

//        static castem_hho_status initialize_workspace(
//                double *const workspace, // worksapce
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** INITIALIZING WORKSPACE" << std::endl;
//                // --- BUILDING THE ELEMENT
//                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
//                std::cout << "- VERICES AS EXTRACTED :" << std::endl;
//                std::cout << verts << std::endl;
//                std::cout << "- CONNECTIVITY AS EXTRACTED :" << std::endl;
//                std::cout << ordening << std::endl;
//                Cast3MElement elem = Cast3MElement(verts, ordening);
//                // --- FILLING THE ELEMENT WORKSPACE FOR GRADIENTS
//                const int offset_gradients = 1;
//                for (int i = 0; i < nq; ++i) {
//                    for (int j = 0; j < grad_size; ++j) {
//                        for (int kk = 0; kk < elem_size; ++kk) {
//                            int array_index = offset_gradients + i * grad_size * elem_size + j * elem_size + kk;
//                            workspace[array_index] += elem.gradient_operators[i](j, kk);
//                        }
//                    }
//                }
//                // --- FILLING THE ELEMENT WORKSPACE FOR THE STABILIZATION
//                const int offset_stabilization = offset_gradients + nq * grad_size * elem_size;
//                for (int i = 0; i < elem_size; ++i) {
//                    for (int j = 0; j < elem_size; ++j) {
//                        int array_index = offset_stabilization + i * elem_size + j;
//                        workspace[array_index] += elem.stabilization_operator(i, j);
//                    }
//                }
//                std::cout << "******** WORKSPACE INITIALIZED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }
//
//        static std::array<EigMat<grad_size, elem_size>, nq> retrieve_gradient_operators(
//                double *const workspace // worksapce
//        ) {
//            std::cout << "--- RETRIEVING GRADIENT OPERATORS" << std::endl;
//            // --- BUILDING THE ELEMENT
//            std::array<EigMat<grad_size, elem_size>, nq> grads;
//            // --- FILLING THE ELEMENT WORKSPACE FOR GRADIENTS
//            const int offset_gradients = 1;
//            for (int i = 0; i < nq; ++i) {
//                for (int j = 0; j < grad_size; ++j) {
//                    for (int kk = 0; kk < elem_size; ++kk) {
//                        int array_index = offset_gradients + i * grad_size * elem_size + j * elem_size + kk;
//                        grads[i](j, kk) = workspace[array_index];
//                    }
//                }
//                std::cout << "- GRADIENT OPERATOR FROM WORKSPACE AT GAUSS POINT " << i << " :" << std::endl;
//                grads[i].print();
//            }
//            std::cout << "--- GRADIENT OPERATORS RETRIEVED" << std::endl;
//            return grads;
//        }
//
//        static EigMat<elem_size, elem_size> retrieve_stabilization(
//                double *const workspace // worksapce
//        ) {
//            std::cout << "--- RETRIEVING STABILIZATION OPERATOR" << std::endl;
//            EigMat<elem_size, elem_size> stab;
//            const int offset_gradients = 1;
//            const int offset_stabilization = offset_gradients + nq * grad_size * elem_size;
//            for (int i = 0; i < elem_size; ++i) {
//                for (int j = 0; j < elem_size; ++j) {
//                    int array_index = offset_stabilization + i * elem_size + j;
//                    stab(i, j) = workspace[array_index];
//                }
//            }
//            std::cout << "- STABILIZATION OPERATOR FROM WORKSPACE :" << std::endl;
//            stab.print();
//            std::cout << "--- STABILIZATION OPERATOR RETRIEVED" << std::endl;
//            return stab;
//        }
//
//        static castem_hho_status get_gauss_data(
//                double *const pts, // gauss points
//                double *const wts, // gauss weights
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** WRITING GAUSS DATA" << std::endl;
//                const EigMap<2, Cast3MElement::num_nodes> verts = EigMap<2, Cast3MElement::num_nodes>(gd->vertices_coordinates);
//                const EigMapIntR<2, Cast3MElement::num_nodes> ordening = EigMapIntR<2, Cast3MElement::num_nodes>(gd->connectivity);
//                Cast3MElement elem = Cast3MElement(verts, ordening);
//                std::cout << "--- FILLING GAUSS OUTPUT" << std::endl;
//                for (int i = 0; i < nq; ++i) {
//                    EigMat<2, 1> qp = elem.cell_cmpt.getQuadPoint(i);
//                    real qw = elem.cell_cmpt.getQuadWeight(i);
//                    int pos_x = i;
//                    int pos_y = nq + i;
//                    pts[pos_x] = qp(0);
//                    pts[pos_y] = qp(1);
//                    wts[i] = qw;
//                }
//                std::cout << "******** GAUSS DATA WRITTEN" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }
//
//        static castem_hho_status compute_gradients(
//                double *const workspace, //worksapce
//                double *const gradient_data, // OUT
//                const double *const faces_unknows, // IN
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** COMPUTING GRADIENTS" << std::endl;
//                std::array<EigMat<grad_size, elem_size>, nq> grad_ops = retrieve_gradient_operators(workspace);
//                EigMat<face_block_size, 1> faces_u;
//                EigMat<elem_size, 1> elem_u;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    int pos = getWorkSpaceSize() + i;
//                    elem_u(i) = workspace[pos];
//                }
//                for (int i = cell_block_size; i < face_block_size; ++i) {
//                    elem_u(i) = faces_unknows[i];
//                }
//                std::cout << "--- FILLING GRADIENTS OUTPUT" << std::endl;
//                for (int i = 0; i < nq; ++i) {
//                    EigMat<1, grad_size> grad = grad_ops[i] * elem_u;
//                    for (int j = 0; j < grad_size; ++j) {
//                        int pos = grad_size * i + j;
//                        gradient_data[pos] = grad(j);
//                    }
//                }
//                std::cout << "******** GRADIENTS COMPUTED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }
//
//        static castem_hho_status compute_internal_forces(
//                double *const workspace, //worksapce
//                double *const internal_forces_data, // OUT
//                const double *const stress, // IN
//                const double *const faces_unknows, // IN
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** COMPUTING INTERNAL FORCES" << std::endl;
//                std::array<EigMat<grad_size, elem_size>, nq> grad_ops = retrieve_gradient_operators(workspace);
//                EigMat<elem_size, elem_size> stab = retrieve_stabilization(workspace);
//                EigMat<face_block_size, 1> faces_u;
//                EigMat<elem_size, 1> intforce;
////                for (int i = 0; i < face_block_size; ++i) {
////                    faces_u(i) = faces_unknows[i];
////                }
//                EigMat<elem_size, 1> elem_u;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    int pos = getWorkSpaceSize() + i;
//                    elem_u(i) = workspace[pos];
//                }
//                for (int i = cell_block_size; i < face_block_size; ++i) {
//                    elem_u(i) = faces_unknows[i];
//                }
//                for (int i = 0; i < nq; ++i) {
//                    EigMat<grad_size, 1> stress_vec;
//                    for (int j = 0; j < grad_size; ++j) {
//                        int pos = grad_size * i + j;
//                        stress_vec(j) = stress[pos];
//                    }
//                    intforce += grad_ops[i].transpose() * stress_vec + stab * elem_u;
//                }
//                std::cout << "--- FILLING INTERNAL FORCES OUTPUT" << std::endl;
//                for (int i = 0; i < elem_size; ++i) {
//                    internal_forces_data[i] = intforce(i);
//                }
//                std::cout << "******** INTERNAL FORCES COMPUTED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }
//
//        static castem_hho_status compute_system(
//                double *const workspace, //worksapce
//                double *const system_data, // OUT
//                const double *const ktan, // IN
//                const double *const res, // IN
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** COMPUTING SYSTEM" << std::endl;
//                std::array<EigMat<grad_size, elem_size>, nq> grad_ops = retrieve_gradient_operators(workspace);
//                EigMat<elem_size, elem_size> stab = retrieve_stabilization(workspace);
//                EigMat<elem_size, 1> residual;
//                for (int i = 0; i < face_block_size; ++i) {
//                    residual(i) = res[i];
//                }
//                EigMat<elem_size, elem_size> Kelem = EigMat<elem_size, elem_size>::Zero();
//                for (int i = 0; i < nq; ++i) {
//                    EigMat<grad_size, grad_size> tangent_op;
//                    for (int j = 0; j < grad_size; ++j) {
//                        for (int kk = 0; kk < grad_size; ++kk) {
//                            int pos = j * grad_size + kk;
//                            tangent_op(j, kk) = ktan[pos];
//                        }
//                    }
//                    std::cout << "!!! SETTING IDENTITY TO TANGENT OPERATOR FOR INVERTIBILITY --> TO BE REMOVED WITH ACTUAL VALUES" << std::endl;
//                    tangent_op = EigMat<grad_size, grad_size>::Identity();
//                    Kelem += grad_ops[i].transpose() * tangent_op * grad_ops[i] + stab;
//                }
//                EigMat<cell_block_size, cell_block_size> m_cell_cell = Kelem.template block<cell_block_size, cell_block_size>(0,0);
//                EigMat<cell_block_size, face_block_size> m_cell_faces = Kelem.template block<cell_block_size, face_block_size>(0,cell_block_size);
//                EigMat<face_block_size, cell_block_size> m_faces_cell = Kelem.template block<face_block_size, cell_block_size>(cell_block_size,0);
//                EigMat<face_block_size, face_block_size> m_faces_faces = Kelem.template block<face_block_size, face_block_size>(cell_block_size,cell_block_size);
//                EigMat<face_block_size, 1> v_faces = residual.template tail<face_block_size>();
//                EigMat<cell_block_size, 1> v_cell = residual.template head<cell_block_size>();
//                std::cout << "--- INVERTING CELL BLOCK MATRIX" << std::endl;
//                EigMat<cell_block_size, cell_block_size> cell_eye = EigMat<cell_block_size, cell_block_size>::Identity();
//                EigMat<cell_block_size, cell_block_size> m_cell_cell_inv = m_cell_cell.llt().template solve(cell_eye);
//                std::cout << "--- CONDENSATING" << std::endl;
//                EigMat<face_block_size, face_block_size> K_cond = m_faces_faces - ((m_faces_cell * m_cell_cell_inv) * m_cell_faces);
//                EigMat<face_block_size, 1> R_cond = v_faces - (m_faces_cell * m_cell_cell_inv) * v_cell;
//                std::cout << "--- FILLING K_COND OUTPUT (COLUMN WISE)" << std::endl;
//                for (int i = 0; i < face_block_size; ++i) {
//                    for (int j = 0; j < face_block_size; ++j) {
//                        int pos = i * face_block_size + j;
//                        system_data[pos] = K_cond(i,j);
//                    }
//                }
//                std::cout << "--- FILLING R_COND OUTPUT" << std::endl;
//                for (int i = 0; i < face_block_size; ++i) {
//                    int pos = face_block_size * face_block_size + i;
//                    system_data[pos] = R_cond(i);
//                }
//                std::cout << "--- UPDATING CELL BLOCK INVERT MATRIX IN WORKSPACE" << std::endl;
//                const int offset_stabilization = 1 + nq * grad_size * elem_size + elem_size * elem_size;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    for (int j = 0; j < cell_block_size; ++j) {
//                        int pos = offset_stabilization + i * cell_block_size + j;
//                        workspace[pos] = m_cell_cell_inv(i,j);
//                    }
//                }
//                std::cout << "--- UPDATING CELL FACE COUPLING BLOCK MATRIX IN WORKSPACE" << std::endl;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    for (int j = 0; j < face_block_size; ++j) {
//                        int pos = offset_stabilization + cell_block_size * cell_block_size + i * face_block_size + j;
//                        workspace[pos] = m_cell_faces(i,j);
//                    }
//                }
//                std::cout << "--- UPDATING CELL RESIDUAL IN WORKSPACE" << std::endl;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    int pos = offset_stabilization + cell_block_size * cell_block_size + face_block_size * cell_block_size;
//                    workspace[pos] = v_cell(i);
//                }
//                std::cout << "******** SYSTEM COMPUTED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }
//
//        static EigMat<cell_block_size, cell_block_size> retrieve_m_cell_cell_inv(
//                double *const workspace // worksapce
//        ) {
//            std::cout << "--- RETRIEVING CELL BLOCK INVERT MATRIX" << std::endl;
//            EigMat<cell_block_size, cell_block_size> m_cell_cell_inv;
//            const int offset_stabilization = 1 + nq * grad_size * elem_size + elem_size * elem_size;
//            for (int i = 0; i < cell_block_size; ++i) {
//                for (int j = 0; j < cell_block_size; ++j) {
//                    int array_index = offset_stabilization + i * cell_block_size + j;
//                    m_cell_cell_inv(i, j) = workspace[array_index];
//                }
//            }
//            std::cout << "- CELL BLOCK INVERT MATRIX FROM WORKSPACE :" << std::endl;
//            m_cell_cell_inv.print();
//            std::cout << "--- CELL BLOCK INVERT MATRIX RETRIEVED" << std::endl;
//            return m_cell_cell_inv;
//        }
//
//        static EigMat<cell_block_size, face_block_size> retrieve_m_cell_faces(
//                double *const workspace // worksapce
//        ) {
//            std::cout << "--- RETRIEVING CELL FACE COUPLING BLOCK MATRIX" << std::endl;
//            EigMat<cell_block_size, face_block_size> m_cell_faces;
//            const int offset_stabilization = 1 + nq * grad_size * elem_size + elem_size * elem_size + cell_block_size * cell_block_size;
//            for (int i = 0; i < cell_block_size; ++i) {
//                for (int j = 0; j < face_block_size; ++j) {
//                    int pos = offset_stabilization + i * face_block_size + j;
//                    m_cell_faces(i,j) = workspace[pos];
//                }
//            }
//            std::cout << "- CELL FACE COUPLING BLOCK MATRIX FROM WORKSPACE :" << std::endl;
//            m_cell_faces.print();
//            std::cout << "--- CELL FACE COUPLING BLOCK MATRIX RETRIEVED" << std::endl;
//            return m_cell_faces;
//        }
//
//        static EigMat<cell_block_size, 1> retrieve_cell_residual(
//                double *const workspace // worksapce
//        ) {
//            std::cout << "--- RETRIEVING CELL RESIDUAL" << std::endl;
//            EigMat<cell_block_size, 1> cell_residual;
//            const int offset_stabilization = 1 + nq * grad_size * elem_size + elem_size * elem_size + cell_block_size * cell_block_size + cell_block_size * face_block_size;
//            for (int i = 0; i < cell_block_size; ++i) {
//                int pos = offset_stabilization + i;
//                cell_residual(i) = workspace[pos];
//            }
//            std::cout << "- CELL RESIDUAL FROM WORKSPACE :" << std::endl;
//            cell_residual.print();
//            std::cout << "--- CELL RESIDUAL RETRIEVED" << std::endl;
//            return cell_residual;
//        }
//
//        static castem_hho_status update_cell_unknows(
//                double *const workspace, //worksapce
//                const double *const faces_correction, // IN
//                const castem_hho_element_geometry *const gd
//        ) {
//            try {
//                std::cout << "******** UPDATING CELL UNKNOWNS" << std::endl;
//                EigMat<cell_block_size, cell_block_size> mccinv = retrieve_m_cell_cell_inv(workspace);
//                EigMat<cell_block_size, face_block_size> mcf = retrieve_m_cell_faces(workspace);
//                EigMat<cell_block_size, 1> vc = retrieve_cell_residual(workspace);
//                EigMat<face_block_size, 1> face_corr;
//                for (int i = 0; i < face_block_size; ++i) {
//                    face_corr(i) = faces_correction[i];
//                }
//                EigMat<cell_block_size, 1> cell_corr = mccinv * (vc - mcf * face_corr);
//                const int offset_stabilization = 1 + nq * grad_size * elem_size + elem_size * elem_size + cell_block_size * cell_block_size + cell_block_size * face_block_size + cell_block_size;
//                for (int i = 0; i < cell_block_size; ++i) {
//                    int pos = offset_stabilization + i;
//                    workspace[pos] += cell_corr(i);
//                }
//                std::cout << "******** CELL UNKNOWNS UPDATED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }

//        static castem_hho_status getDescription(
//                castem_hho_element_description *const d,
//                const castem_hho_element_geometry *const gd
//        ){
//            try {
//                std::cout << "******** INITIALIZING ELEMENT DESCRIPTION" << std::endl;
//                d->element_workspace_size = getWorkSpaceSize();
//                d->get_gauss_data = get_gauss_data;
//                d->initialize_workspace = initialize_workspace;
//                d->compute_gradients = compute_gradients;
//                d->compute_internal_forces = compute_internal_forces;
//                d->compute_system = compute_system;
//                d->decondensate = update_cell_unknows;
//                std::cout << "******** ELEMENT DESCRIPTION INITIALIZED" << std::endl;
//            } catch (...) {
//                return castem_hho_handle_cxx_exception();
//            }
//            return castem_hho_report_success();
//        }

    };

}

#endif //FETA_FUNCTIONS_HXX
