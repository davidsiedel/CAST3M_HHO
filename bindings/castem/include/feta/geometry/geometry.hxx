//
// Created by dsiedel on 24/06/2021.
//

#ifndef FETA_GEOMETRY_HXX
#define FETA_GEOMETRY_HXX

#include <utility>

#include "feta/feta.hxx"

namespace feta::geometry{

    struct Node{
        Node(intg tag, const EigMat<3, 1> &vector, intg cell_weight, intg face_weight) :
                tag(tag),
                coordinates(vector),
                cell_weight(cell_weight),
                face_weight(face_weight) {}

        const intg tag;
        const EigMat<3, 1> coordinates;
        const intg cell_weight;
        const intg face_weight;
    };

    struct QuadraturePoint{
        QuadraturePoint(intg tag, const EigMat<3, 1> &vector, real weight) :
                tag(tag),
                coordinates(vector),
                weight(weight) {}

        template<intg dim>
        static QuadraturePoint fromVector(intg tag, const EigMat<dim, 1> &vector, real weight) {
            EigMat<3, 1> coordinates = EigMat<3, 1>::Zero();
            for (int i = 0; i < dim; ++i) {
                coordinates(i) = vector(i);
            }
            return QuadraturePoint(tag, coordinates, weight);
        }

        static QuadraturePoint fromCoordinates(intg tag, real x, real weight) {
            EigMat<3, 1> coordinates = EigMat<3, 1>::Zero();
            coordinates(0) = x;
            return QuadraturePoint(tag, coordinates, weight);
        }

        static QuadraturePoint fromCoordinates(intg tag, real x, real y, real weight) {
            EigMat<3, 1> coordinates = EigMat<3, 1>::Zero();
            coordinates(0) = x;
            coordinates(1) = y;
            return QuadraturePoint(tag, coordinates, weight);
        }

        static QuadraturePoint fromCoordinates(intg tag, real x, real y, real z, real weight) {
            EigMat<3, 1> coordinates = EigMat<3, 1>::Zero();
            coordinates(0) = x;
            coordinates(1) = y;
            coordinates(2) = z;
            return QuadraturePoint(tag, coordinates, weight);
        }

        const intg tag;
        const EigMat<3, 1> coordinates;
        const real weight;
    };

    template<intg d, intg n>
    struct Nodes{

        explicit Nodes(EigMat<d, n> coordinates) : coordinates(coordinates) {}

        static constexpr intg dim = d;
        static constexpr intg num = n;
        const EigMat<d, n> coordinates;

        static Nodes<d, n> fromArray(real *a) {

        }

        EigMat<d, d> getRotationMatrix() const {
            if constexpr(d == 2) {
//            if (d == 2) {
                return getCurveRotationMatrix();
            } else if constexpr(d == 3) {
//            } else if (d == 3) {
                return getSurfaceRotationMatrix();
            }
        }

        EigMat<d - 1, n> getHyperplanarCoordinates() const {
            return (this->getRotationMatrix() * coordinates).template block<d - 1, n>(0, 0);
        }

        Nodes<d - 1, n> getHyperplanarNodes() const {
            return Nodes<d - 1, n>(getHyperplanarCoordinates());
        }

        EigMat<d, 1> getBarycenter() const {
            EigMat<d, 1> barycenter = EigMat<d, 1>::Zero();
            for (int i = 0; i < n; ++i) {
                barycenter += coordinates.col(i);
            }
            barycenter /= double(n);
            return barycenter;
        }

        EigMat<d, 1> getBounds() const {
            EigMat<d, 1> bounds = EigMat<d, 1>::Zero();
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (j > i) {
                        for (int k = 0; k < d; ++k) {
                            real quant = abs(coordinates(k, i) - coordinates(k, j));
                            if (quant > abs(bounds(k))) {
                                bounds(k) = quant;
                            }
                        }
                    }
                }
            }
            return bounds;
        }

    private:

        EigMat<2, 2> getCurveRotationMatrix() const {
            EigMat<2, 2> rotation_matrix;
            EigMat<2, 1> edge = coordinates.col(1) - coordinates.col(0);
            edge /= edge.norm();
            rotation_matrix(0, 0) = edge(0);
            rotation_matrix(0, 1) = edge(1);
            rotation_matrix(1, 0) = edge(1);
            rotation_matrix(1, 1) = -edge(0);
            return rotation_matrix;
        }

        EigMat<3, 3> getSurfaceRotationMatrix() const {
            EigMat<3, 3> rotation_matrix;
            EigMat<3, 1> edge_0 = coordinates.col(2) - coordinates.col(0);
            EigMat<3, 1> edge_t = coordinates.col(1) - coordinates.col(0);
            EigMat<3, 1> edge_2 = edge_0.cross(edge_t);
            EigMat<3, 1> edge_1 = edge_2.cross(edge_0);
            edge_0 = edge_0 / edge_0.norm();
            edge_1 = edge_1 / edge_1.norm();
            edge_2 = edge_2 / edge_2.norm();
            rotation_matrix.row(0) = edge_0;
            rotation_matrix.row(1) = edge_1;
            rotation_matrix.row(2) = edge_2;
            return rotation_matrix;
        }
    };

    template<intg dim_eucli, intg n_nodes>
    void checkNodesConsistency(){
        if constexpr(n_nodes < dim_eucli) throw std::logic_error("unsufficient number fo nodes");
//        if (n_nodes < dim_eucli) throw std::logic_error("unsufficient number fo nodes");
    }

//    template<intg dim_eucli>
//    EigMat<3, 1> getVector(const EigMat<dim_eucli, 1> v) {
//        EigMat<3, 1> qp = EigMat<3, 1>::Zero();
//        for (int i = 0; i < dim_eucli; ++i) {
//            qp(i) = v(i);
//        }
//        return qp;
//    }

    template<intg n_nodes>
    EigMat<2, 2> getRotationMatrix(const EigMat<2, n_nodes> &nodes){
        checkNodesConsistency<2, n_nodes>();
        EigMat<2, 2> rotation_matrix;
        EigMat<2, 1> edge = nodes.col(1) - nodes.col(0);
        edge /= edge.norm();
        rotation_matrix(0, 0) = edge(0);
        rotation_matrix(0, 1) = edge(1);
        rotation_matrix(1, 0) = edge(1);
        rotation_matrix(1, 1) = -edge(0);
        return rotation_matrix;
    }

    template<intg n_nodes>
    EigMat<3, 3> getRotationMatrix(const EigMat<3, n_nodes> &nodes){
        checkNodesConsistency<3, n_nodes>();
        EigMat<3, 3> rotation_matrix;
        EigMat<3, 1> edge_0 = nodes.col(2) - nodes.col(0);
        EigMat<3, 1> edge_t = nodes.col(1) - nodes.col(0);
        EigMat<3, 1> edge_2 = edge_0.cross(edge_t);
        EigMat<3, 1> edge_1 = edge_2.cross(edge_0);
        edge_0 = edge_0 / edge_0.norm();
        edge_1 = edge_1 / edge_1.norm();
        edge_2 = edge_2 / edge_2.norm();
        rotation_matrix.row(0) = edge_0;
        rotation_matrix.row(1) = edge_1;
        rotation_matrix.row(2) = edge_2;
        return rotation_matrix;
    }

    template<intg dim_euclidean, intg n_nodes>
    EigMat<dim_euclidean, 1> getBarycenter(const EigMat<dim_euclidean, n_nodes> &nodes) {
        EigMat<dim_euclidean, 1> barycenter = EigMat<dim_euclidean, 1>::Zero();
        for (int i = 0; i < n_nodes; ++i) {
            barycenter += nodes.col(i);
        }
        barycenter /= double(n_nodes);
        return barycenter;
    }

    template<intg dim_euclidean, intg n_nodes>
    EigMat<dim_euclidean, 1> getBounds(const EigMat<dim_euclidean, n_nodes> &nodes) {
        EigMat<dim_euclidean, 1> bounds = EigMat<dim_euclidean, 1>::Zero();
        for (int i = 0; i < n_nodes; ++i) {
            for (int j = 0; j < n_nodes; ++j) {
                if (j > i) {
                    for (int k = 0; k < dim_euclidean; ++k) {
                        real quant = abs(nodes(k, i) - nodes(k, j));
                        if (quant > abs(bounds(k))) {
                            bounds(k) = quant;
                        }
                    }
                }
            }
        }
        return bounds;
    }

//    template<intg dim_eucli, intg n_nodes>
//    EigMat<dim_eucli - 1, n_nodes> getNodesInHyperplane(const EigMat<dim_eucli, n_nodes> &nodes) {
//        const EigMat<dim_eucli, dim_eucli> rotation_matrix = getRotationMatrix2<dim_eucli, n_nodes>(nodes);
//        return (rotation_matrix * nodes).template block<dim_eucli - 1, n_nodes>(0, 0);
//    }
}

#endif //FETA_GEOMETRY_HXX
