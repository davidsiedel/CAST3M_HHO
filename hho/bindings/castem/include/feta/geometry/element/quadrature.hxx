//
// Created by dsiedel on 25/06/2021.
//

#ifndef FETA_QUADRATURE_HXX
#define FETA_QUADRATURE_HXX

#include "feta/geometry/element/element.hxx"

namespace feta::geometry {

    template<intg d, typename Element_t>
    struct StandardQuadrature {
        explicit StandardQuadrature(const Nodes<d, Element_t::num_nodes> &nodes) :
                coordinates(get(nodes)),
                weights(g(nodes)) {}

        std::vector<QuadraturePoint> getStandardElementQuadraturePoints() const {
            std::vector<QuadraturePoint> quada;
            quada.reserve(dim);
            for (int i = 0; i < dim; ++i) {
                QuadraturePoint qpa = QuadraturePoint::fromVector<d>(0, coordinates.col(i), weights(i));
                quada.push_back(qpa);
            }
            return quada;
        }

    private:

        using Table = typename Element_t::Table;
        static constexpr intg dim = Table::dim_quadrature;

        static constexpr intg dim_e = Element_t::dim_element;
        static constexpr intg num_nodes = Element_t::num_nodes;


        const EigMat<d, dim> coordinates;
        const EigMat<1, dim> weights;

        EigMat<1, dim> g(const Nodes<dim_e, num_nodes> &nodes) const {
            static constexpr intg dq = Table::dim_quadrature;
            static constexpr intg dim_jac = dim_e * dim_e;
            EigMat<1, dq> quadrature_weights;
            const auto quadrature_reference_weights = EigMapC<dq, 1>(Table::weights.data());
            for (int i = 0; i < dq; ++i) {
                const real reference_weight = quadrature_reference_weights(i);
                const auto jacobian_matrix = EigMap<dim_jac, num_nodes>(Table::jacobians[i].data());
                EigMat<dim_e, dim_e> jacc;
                for (int j = 0; j < dim_e; ++j) {
                    for (int k = 0; k < dim_e; ++k) {
                        int a = dim_e * j + k;
                        jacc(j, k) = jacobian_matrix.row(a).dot(nodes.coordinates.row(j));
                    }
                }
                const real weight = reference_weight * abs(jacc.determinant());
                quadrature_weights(i) = weight;
            }
            return quadrature_weights;
        }

        EigMat<1, dim> g(const Nodes<dim_e + 1, num_nodes> &nodes) const {
            Nodes<dim_e, num_nodes> vcs = nodes.getHyperplanarNodes();
            return g(vcs);
        }

        EigMat<d, dim> get(const Nodes<d, num_nodes> &nodes) const {
            EigMap<dim, num_nodes> quadrature_reference_points = EigMap<dim, num_nodes>(Table::points.data());
            EigMat<d, num_nodes> nds = nodes.coordinates;
            EigMat<d, dim> quadrature_points = (quadrature_reference_points * nds.transpose()).transpose();
            return quadrature_points;
        }
    };

    template<intg d, ElementShape elem_shp, intg num_nodes, intg k>
    struct GenericQuadrature {
        explicit GenericQuadrature(const Nodes<d, num_nodes> &nodes) : qps(setQuadraturePoints(nodes)) {}

        std::vector<QuadraturePoint> getQuadraturePoints() {
            return qps;
        }
        using Elem_t = Element<elem_shp, num_nodes, k>;
        static constexpr intg dim = Elem_t::dim_quadrature;
    private:
        const std::vector<QuadraturePoint> qps;

        std::vector<QuadraturePoint> setQuadraturePoints(const Nodes<d, num_nodes> &nodes) {
            return StandardQuadrature<d, Elem_t>(nodes).getStandardElementQuadraturePoints();
        }
    };

    template<intg d, intg num_nodes, intg k>
    struct GenericQuadrature<d, ElementShape::Polygon, num_nodes, k> {
        explicit GenericQuadrature(const Nodes<d, num_nodes> &nodes) : qps(setQuadraturePoints(nodes)) {}

        std::vector<QuadraturePoint> getQuadraturePoints() {
            return qps;
        }
        using Elem_t = Element<ElementShape::Polygon, num_nodes, k>;
        static constexpr intg dim = Elem_t::dim_quadrature;

    private:

        const std::vector<QuadraturePoint> qps;

        using Tri_t = Element<ElementShape::Triangle, 3, k>;
        static constexpr intg dqt = Elem_t::Table::dim_quadrature;

        static std::vector<QuadraturePoint> setQuadraturePoints(const Nodes<d, num_nodes> &nodes) {
            std::array<EigMat<d, 3>, num_nodes> partition = getPolygonPartition(nodes);
            std::vector<QuadraturePoint> quada;
            quada.reserve(dim);
            for (int i = 0; i < num_nodes; ++i) {
                Nodes<d, 3> tri = Nodes<d, 3>(partition[i]);
                std::vector<QuadraturePoint> quads = StandardQuadrature<d, Tri_t>(tri).getStandardElementQuadraturePoints();
                // -- moving values
                for (QuadraturePoint qpp: quads) {
                    quada.push_back(qpp);
                }
//                quada.insert(quada.end(), std::make_move_iterator(quads.begin()), std::make_move_iterator(quads.end()));
            }
            return quada;
        }

        static std::array<EigMat<d, 3>, num_nodes> getPolygonPartition(const Nodes<d, num_nodes> &nodes) {
            const EigMat<d, 1> barycenter = nodes.getBarycenter();
//            std::array<Nodes<d, 3>, num_nodes> partition;
            std::array<EigMat<d, 3>, num_nodes> partition;
            for (int i = 0; i < num_nodes; ++i) {
                EigMat<d, 3> triangle_nodes;
                if (i != num_nodes - 1) {
                    triangle_nodes.col(0) = nodes.coordinates.col(i);
                    triangle_nodes.col(2) = nodes.coordinates.col(i + 1);
                    triangle_nodes.col(1) = barycenter;
                } else {
                    triangle_nodes.col(0) = nodes.coordinates.col(i);
                    triangle_nodes.col(2) = nodes.coordinates.col(0);
                    triangle_nodes.col(1) = barycenter;
                }
//                partition[i] = Nodes<d, 3>(triangle_nodes);
                partition[i] = triangle_nodes;
            }
            return partition;
        }

    };

//    template<intg d, ElementType elem_t, intg k>
//    struct Quadr2 : public GenericQuadrature<d, Elem2<elem_t, k>::shape, Elem2<elem_t, k>::num_nodes, k> {
//        explicit Quadr2(const Nodes<d, Elem2<elem_t, k>::num_nodes> &nodes) : GenericQuadrature<d, Elem2<elem_t, k>::shape, Elem2<elem_t, k>::num_nodes, k>(nodes) {}
//    };
//    template<intg d, ElementType elem_t, intg k>
//    std::vector<QuadraturePoint> getQuadraturePoints(const Nodes<d, Elem2<elem_t, k>::num_nodes> &nodes) {
//        return GenericQuadrature<d, Elem2<elem_t, k>::shape, Elem2<elem_t, k>::num_nodes, k>(nodes).getQuadraturePoints();
//    }

    template<intg d, ElementType elem_t, intg k>
    struct Quadrature {
        explicit Quadrature(const Nodes<d, Elem2<elem_t, k>::num_nodes> &nodes) : qps(setQuadraturePoints(nodes)) {}

        std::vector<QuadraturePoint> getQuadraturePoints() const {
            return qps;
        }
        using Quada = GenericQuadrature<d, Elem2<elem_t, k>::shape, Elem2<elem_t, k>::num_nodes, k>;
        static constexpr intg dim = Quada::dim;

    private:

        const std::vector<QuadraturePoint> qps;

        std::vector<QuadraturePoint> setQuadraturePoints(const Nodes<d, Elem2<elem_t, k>::num_nodes> &nodes) {
            return Quada(nodes).getQuadraturePoints();
        };
    };
}

#endif //FETA_QUADRATURE_HXX
