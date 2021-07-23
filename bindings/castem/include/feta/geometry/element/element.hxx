//
// Created by dsiedel on 25/06/2021.
//

#ifndef FETA_ELEMENT_HXX
#define FETA_ELEMENT_HXX

#include "feta/geometry/geometry.hxx"
#include "feta/geometry/element/gauss_quadrature/gauss_linear_segment.hxx"
#include "feta/geometry/element/gauss_quadrature/gauss_linear_triangle.hxx"
#include "feta/geometry/element/gauss_quadrature/gauss_linear_quadrangle.hxx"

namespace feta::geometry{

    enum struct ElementShape {
        Segment,
        Triangle,
        Quadrangle,
        Polygon
    };

    template<ElementShape elem_t, intg n_nodes, intg k>
    struct Element;

    template<intg k>
    struct Element<ElementShape::Segment, 2, k>{
        using Table = GaussLinearSegment<k>;
        using Reference = LinearSegment;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = Reference::dim_element;
        static constexpr intg dim_jacobian = Reference::dim_jacobian;
        static constexpr intg num_nodes = Reference::num_nodes;
    };

    template<intg k>
    struct Element<ElementShape::Triangle, 3, k>{
        using Table = GaussLinearTriangle<k>;
        using Reference = LinearTriangle;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = Reference::dim_element;
        static constexpr intg dim_jacobian = Reference::dim_jacobian;
        static constexpr intg num_nodes = Reference::num_nodes;
    };

    template<intg k>
    struct Element<ElementShape::Quadrangle, 4, k>{
        using Table = GaussLinearQuadrangle<k>;
        using Reference = LinearQuadrangle;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = Reference::dim_element;
        static constexpr intg dim_jacobian = Reference::dim_jacobian;
        static constexpr intg num_nodes = Reference::num_nodes;
    };

    template<intg n_nodes, intg k>
    struct Element<ElementShape::Polygon, n_nodes, k>{
        using Table = GaussLinearTriangle<k>;
        using Reference = LinearTriangle;
        static constexpr intg dim_quadrature = Table::dim_quadrature * n_nodes;
        static constexpr intg dim_element = Reference::dim_element;
        static constexpr intg dim_jacobian = Reference::dim_jacobian;
        static constexpr intg num_nodes = n_nodes;
    };

    enum struct ElementType {
        LinearSegment,
        LinearTriangle,
        LinearQuadrangle,
        LinearPentagon,
        LinearHexagon,
        LinearHeptagon,
        LinearOctagon,
        LinearEnneagon,
        LinearDecagon,
        LinearHendecagon,
        LinearDodecagon,
        LinearTridecagon,
        LinearTetradecagon,
        LinearPentadecagon,
        LinearHexadecagon,
        LinearHeptadecagon,
        LinearOctadecagon,
        LinearEnneadecagon,
        LinearIcosagon
    };

    template<ElementType elem_t, intg k>
    struct Elem2;

    template<intg  k>
    struct Elem2<ElementType::LinearSegment, k>{
        using Table = GaussLinearSegment<k>;
        static constexpr ElementShape shape = ElementShape::Segment;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = 1;
        static constexpr intg dim_jacobian = 1;
        static constexpr intg num_nodes = 2;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearTriangle, k>{
        using Table = GaussLinearTriangle<k>;
        static constexpr ElementShape shape = ElementShape::Triangle;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = 2;
        static constexpr intg dim_jacobian = 4;
        static constexpr intg num_nodes = 3;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearQuadrangle, k>{
        using Table = GaussLinearQuadrangle<k>;
        static constexpr ElementShape shape = ElementShape::Quadrangle;
        static constexpr intg dim_quadrature = Table::dim_quadrature;
        static constexpr intg dim_element = 2;
        static constexpr intg dim_jacobian = 4;
        static constexpr intg num_nodes = 4;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearPentagon, k>{
        using Table = GaussLinearTriangle<k>;
        static constexpr ElementShape shape = ElementShape::Polygon;
        static constexpr intg dim_quadrature = Table::dim_quadrature * 5;
        static constexpr intg dim_element = 2;
        static constexpr intg dim_jacobian = 4;
        static constexpr intg num_nodes = 5;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearHexagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 6;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 6;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearHeptagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 7;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 7;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearOctagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 8;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 8;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearEnneagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 9;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 9;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearDecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 10;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 10;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearHendecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 11;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 11;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearDodecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 12;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 12;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearTridecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 13;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 13;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearTetradecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 14;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 14;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearPentadecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 15;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 15;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearHexadecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 16;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 16;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearHeptadecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 17;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 17;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearOctadecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 18;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 18;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearEnneadecagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 19;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 19;
    };

    template<intg  k>
    struct Elem2<ElementType::LinearIcosagon, k>{
      using Table = GaussLinearTriangle<k>;
      static constexpr ElementShape shape = ElementShape::Polygon;
      static constexpr intg dim_quadrature = Table::dim_quadrature * 20;
      static constexpr intg dim_element = 2;
      static constexpr intg dim_jacobian = 4;
      static constexpr intg num_nodes = 20;
    };

    template<typename elem_t, intg k>
    struct ElementData : Element<elem_t::shape, elem_t::num_nodes, k>{

    };
}

#endif //FETA_ELEMENT_HXX
