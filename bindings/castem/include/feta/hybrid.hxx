//
// Created by dsiedel on 24/06/2021.
//

#ifndef FETA_HYBRID_HXX
#define FETA_HYBRID_HXX

#include "feta/geometry/element/quadrature.hxx"
#include "feta/basis/basis.hxx"

using namespace feta::geometry;
using namespace feta::basis;

namespace feta::hybrid{

    enum struct FieldType {
        VECTOR_SMALL_STRAIN,
        VECTOR_SMALL_STRAIN_CAST3M,
        VECTOR_FINITE_STRAIN_CAST3M,
        SCALAR
    };

    enum struct DerivationType {
        FULL,
        SYMMETRIC
    };

    enum struct HHOOrder {
        LOW,
        EQUAL,
        HIGH
    };

//    template<intg indx, intg d0, intg d1, real cof>
//    struct VoigtData{
//        static constexpr intg index = indx;
//        static constexpr intg di = d0;
//        static constexpr intg dj = d1;
//        static constexpr real coef = cof;
//    };
//
//    template<typename index_t, typename di_t, typename dj_t, typename coef_t>
//    struct VoigtData2{
//        static constexpr index_t index;
//        static constexpr di_t di;
//        static constexpr dj_t dj;
//        static constexpr coef_t coef;
//
//        index_t getIndex() {
//            return index;
//        }
//    };

    template<intg d, intg k, intg l>
//    template<intg d, intg l, intg k>
    struct HHOSpace{
        HHOSpace() :
                grad_basis(Monomial<d, k>()),
                cell_basis(Monomial<d, l>()),
                reco_basis(Monomial<d, k+1>()),
                face_basis(Monomial<d-1, k>())
                {}
        Monomial<d, k> grad_basis;
        Monomial<d, l> cell_basis;
        Monomial<d, k+1> reco_basis;
        Monomial<d-1, k> face_basis;
        static constexpr intg ck = Monomial<d, k>::dim;
        static constexpr intg cl = Monomial<d, l>::dim;
        static constexpr intg cr = Monomial<d, k+1>::dim;
        static constexpr intg fk = Monomial<d-1, k>::dim;
        static constexpr intg bld_k = 2 * (Monomial<d, k>::order + 1);
        static constexpr intg cpt_k = 2 * Monomial<d, k>::order;
//        static constexpr intg bld_k = 1;
//        static constexpr intg cpt_k = 1;
        static void print() {
            std::cout << "cell basis dimensions : " << cl << std::endl;
            std::cout << "grad basis dimensions : " << ck << std::endl;
            std::cout << "reconstruction basis dimensions : " << cr << std::endl;
            std::cout << "face basis dimensions : " << fk << std::endl;
            std::cout << "build integration order : " << bld_k << std::endl;
            std::cout << "computation integration order : " << cpt_k << std::endl;
        }
    };

    template<intg dim, FieldType field_type>
    struct FieldPolicy;

    template <>
    struct FieldPolicy<2, FieldType::VECTOR_SMALL_STRAIN> {
        static constexpr DerivationType derivation_type = DerivationType::SYMMETRIC;
        static constexpr intg field_size = 2;
        static constexpr intg gradient_size = 4;
        static constexpr intg voigt_size = 3;
        static constexpr std::array<std::tuple<intg, intg, intg, real>, voigt_size> voigt_data = {
                {
                        {0, 0, 0, 1.0},
                        {1, 1, 1, 1.0},
                        {3, 0, 1, SQRT2}
                }
        };
//        static constexpr std::array<VoigtData2<intg, intg, intg, real>, voigt_size> voigt_data = {
//                {
//                        {0, 0, 0, 1.0},
//                        {1, 1, 1, 1.0},
//                        {3, 0, 1, SQRT2}
//                }
//        };
    };

    template <>
    struct FieldPolicy<2, FieldType::VECTOR_SMALL_STRAIN_CAST3M> {
        static constexpr DerivationType derivation_type = DerivationType::SYMMETRIC;
        static constexpr intg field_size = 2;
        static constexpr intg gradient_size = 9;
        static constexpr intg voigt_size = 4;
        static constexpr std::array<std::tuple<intg, intg, intg, real>, voigt_size> voigt_data = {
                {
                        {0, 0, 0, 1.0},
                        {1, 1, 1, 1.0},
                        {3, 0, 1, 1.0},
                        {4, 1, 0, 1.0}
                }
        };
    };

    template <>
    struct FieldPolicy<2, FieldType::VECTOR_FINITE_STRAIN_CAST3M> {
        static constexpr DerivationType derivation_type = DerivationType::FULL;
        static constexpr intg field_size = 2;
        static constexpr intg gradient_size = 9;
        static constexpr intg voigt_size = 4;
        static constexpr std::array<std::tuple<intg, intg, intg, real>, voigt_size> voigt_data = {
                {
                        {0, 0, 0, 1.0},
                        {1, 1, 1, 1.0},
                        {3, 0, 1, 1.0},
                        {4, 1, 0, 1.0}
                }
        };
    };

    template<intg k>
    struct Cast3MFace{
        explicit Cast3MFace(const EigMat<2, 2> &verts) :
        quadrature_points(setQuadrature(verts)),
        nodes(setNodes(verts)){}
        using FaceQuad = Quadrature<2, ElementType::LinearSegment, k>;
        static constexpr intg dim_q = FaceQuad::dim;
        const std::vector<QuadraturePoint> quadrature_points;
        const Nodes<2, 2> nodes;

        EigMat<2, 1> getQuadPoint(intg i) const {
            return quadrature_points[i].coordinates.template head<2>();
//            return quadrature_points[i].coordinates;
        }

        real getQuadWeight(intg i) const {
            return quadrature_points[i].weight;
        }

        EigMat<1, 1> getQuadPointF(intg i) const {
            EigMat<2, 1> pc = getQuadPoint(i);
            EigMat<1, 1> pf = (nodes.getRotationMatrix() * pc).template head<1>();
            return pf;
        }

        EigMat<2, 1> getBarycenter() const {
            return nodes.getBarycenter();
        }

        EigMat<2, 2> getRotationMatrix() const {
            return nodes.getRotationMatrix();
        }

        EigMat<2, 1> getNormalVector() const {
            return nodes.getRotationMatrix().row(1);
        }

        EigMat<2, 1> getBounds() const {
            return nodes.getBounds();
        }

        EigMat<1, 1> getBarycenterF() const {
            return (nodes.getRotationMatrix() * nodes.getBarycenter()).template head<1>();
        }

        EigMat<1, 1> getBoundsF() const {
            Nodes<1, 2> ndsprj = nodes.getHyperplanarNodes();
            return ndsprj.getBounds();
        }

        void print() const {
            std::cout << "-------------------------------" << std::endl;
            std::cout << "--- NODES" << std::endl;
            std::cout << nodes.coordinates << std::endl;
            std::cout << "integration order : " << k << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "dim quadrature : " << quadrature_points.size() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            for (int i = 0; i < quadrature_points.size(); ++i) {
                std::cout << "quadrature point : " << i << " :" << std::endl;
                std::cout << "weight : " << quadrature_points[i].weight << std::endl;
                std::cout << quadrature_points[i].coordinates << std::endl;
            }
            for (int i = 0; i < quadrature_points.size(); ++i) {
                std::cout << "quadrature point in face : " << i << " :" << std::endl;
                std::cout << getQuadPointF(i) << std::endl;
            }
            std::cout << "-------------------------------" << std::endl;
            std::cout << "normal vector : "<< std::endl;
            std::cout << getNormalVector() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "barycenter : "<< std::endl;
            std::cout << getBarycenter() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "rotation matrix : "<< std::endl;
            std::cout << getRotationMatrix() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "bounds : "<< std::endl;
            std::cout << getBounds() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "barycenter in face : "<< std::endl;
            std::cout << getBarycenterF() << std::endl;
            std::cout << "-------------------------------" << std::endl;
            std::cout << "bounds in face : "<< std::endl;
            std::cout << getBoundsF() << std::endl;
        }

    private:
        std::vector<QuadraturePoint> setQuadrature(const EigMat<2, 2> &verts) const {
            Nodes<2, 2> nds = Nodes<2, 2>(verts);
            return FaceQuad(nds).getQuadraturePoints();
        }

        Nodes<2, 2> setNodes(const EigMat<2, 2> &verts) const {
            return Nodes<2, 2>(verts);
        }
    };

    template<ElementType elem_t, intg k>
    struct Cast3MCell{
        Cast3MCell(const EigMat<2, Elem2<elem_t, k>::num_nodes> &verts) :
                quadrature_points(setQuadrature(verts)),
                nodes(setNodes(verts)){}
        using CellQuad = Quadrature<2, elem_t, k>;
        static constexpr intg dim_q = CellQuad::dim;
        static constexpr intg num_nodes = Elem2<elem_t, k>::num_nodes;

        const std::vector<QuadraturePoint> quadrature_points;
        const Nodes<2, num_nodes> nodes;

        EigMat<2, 1> getQuadPoint(intg i) const {
            return quadrature_points[i].coordinates.template head<2>();
//            return quadrature_points[i].coordinates;
        }

        real getQuadWeight(intg i) const {
            return quadrature_points[i].weight;
        }

        EigMat<2, 1> getBarycenter() const {
            return nodes.getBarycenter();
        }

        EigMat<2, 1> getBounds() const {
            return nodes.getBounds();
        }

        void print() const {
            std::cout << "integration order : " << k << std::endl;
            std::cout << "dim quadrature : " << quadrature_points.size() << std::endl;
            for (int i = 0; i < quadrature_points.size(); ++i) {
                std::cout << "quadrature point : " << i << " :" << std::endl;
                std::cout << "weight : " << quadrature_points[i].weight << std::endl;
                std::cout << quadrature_points[i].coordinates << std::endl;
            }
        }

    private:
        std::vector<QuadraturePoint> setQuadrature(const EigMat<2, num_nodes> &verts) const {
            Nodes<2, num_nodes> nds = Nodes<2, num_nodes>(verts);
            return CellQuad(nds).getQuadraturePoints();
        }

        Nodes<2, num_nodes> setNodes(const EigMat<2, num_nodes> &verts) const {
            return Nodes<2, num_nodes>(verts);
        }
    };

    template<ElementType elem_t, FieldType field_t, intg k, intg l>
    struct Cast3MElem{
        Cast3MElem(const EigMat<2, Elem2<elem_t, k>::num_nodes> &vertices,
                   const EigMat2<2, Elem2<elem_t, k>::num_nodes> &conn)
                   :
                   cell(CellBuild(vertices)),
                   faces(setFaces<FaceBuild>(vertices, conn)),
                   cell_cmpt(CellCmpt(vertices)),
                   faces_cmpt(setFaces<FaceCmpt>(vertices, conn)),
                   space(HHOSpace<2, k, l>()),
                   gradient_operators(getGradientOperators()),
                   stabilization_operator(getStabilizationOperator())
//                   recrhs(getReconstructionOperator())
                   {}
        static constexpr intg d = 2;
        using Space = HHOSpace<d, k, l>;
        using Field = FieldPolicy<d, field_t>;
        static constexpr intg bld_k = Space::bld_k;
        static constexpr intg cpt_k = Space::cpt_k;
        using FaceBuild = Cast3MFace<bld_k>;
        using CellBuild = Cast3MCell<elem_t, bld_k>;
        using FaceCmpt = Cast3MFace<cpt_k>;
        using CellCmpt = Cast3MCell<elem_t, cpt_k>;
        //-
        static constexpr intg grad_size = Field::gradient_size;
        static constexpr intg field_size = Field::field_size;
        static constexpr intg voigt_size = Field::voigt_size;
        static constexpr intg num_nodes = CellCmpt::num_nodes;
        static constexpr intg ck = Space::ck;
        static constexpr intg cl = Space::cl;
        static constexpr intg cr = Space::cr;
        static constexpr intg fk = Space::fk;
        static constexpr intg elem_size = field_size * (cl + num_nodes * fk);
        //-
        const Space space;
        //-
        const std::vector<FaceBuild> faces;
        const CellBuild cell;
        //-
        const std::vector<FaceCmpt> faces_cmpt;
        const CellCmpt cell_cmpt;
        const std::vector<EigMat<grad_size, elem_size>> gradient_operators;
        const EigMat<elem_size, elem_size> stabilization_operator;
//        EigMat<cr * field_size, elem_size> recrhs;
    private:

        template<typename FaceType>
        std::vector<FaceType> setFaces(const EigMat<2, num_nodes> &vertices, const EigMat2<2, num_nodes> &conn) const {
            std::vector<FaceType> fcs;
            fcs.reserve(num_nodes);
            for (int i = 0; i < num_nodes; ++i) {
                EigMat2<2, 1> nds_idx = conn.col(i);
                EigMat<2, 2> fc_vtx;
                fc_vtx.col(0) = vertices.col(nds_idx(0));
                fc_vtx.col(1) = vertices.col(nds_idx(1));
                FaceType fc = FaceType(fc_vtx);
                fcs.push_back(fc);
            }
            return fcs;
        }

        EigMat<cr-1, cr-1> getReconstructionLHSInvert() const {
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            EigMat<cr-1, cr-1> eye = EigMat<cr-1, cr-1>::Identity();
            EigMat<cr, cr> lhs = EigMat<cr, cr>::Zero();
            for (int i = 0; i < CellBuild::dim_q; ++i) {
                for (int j = 0; j < Field::field_size; ++j) {
                    const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
                    const real w_q_c = cell.getQuadWeight(i);
                    const EigMat<cr, 1> d_phi_r = space.reco_basis.GetDerivativeVector(x_q_c, x_c, bdc, j);
                    lhs += w_q_c * d_phi_r * d_phi_r.transpose();
                }
            }
            EigMat<cr-1, cr-1> inv = (lhs.template block<cr - 1, cr - 1>(1, 1)).llt().template solve(eye);
            return inv;
        }

        EigMat<cr - 1, elem_size> getReconstructionRHS(intg di) const {
            EigMat<cr, elem_size> rhs = EigMat<cr, elem_size>::Zero();
            // CELL GEOMETRY
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            EigMat<cr, cl> m_stif = EigMat<cr, cl>::Zero();
            for (int i = 0; i < CellBuild::dim_q; ++i) {
                for (int j = 0; j < Field::field_size; ++j) {
                    const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
                    const real w_q_c = cell.getQuadWeight(i);
                    const EigMat<cl, 1> d_phi_l_j = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, j);
                    const EigMat<cr, 1> d_phi_r_j = space.reco_basis.GetDerivativeVector(x_q_c, x_c, bdc, j);
                    m_stif += w_q_c * d_phi_r_j * d_phi_l_j.transpose();
                }
            }
            rhs.template block<cr, cl>(0, di * cl) += m_stif;
            for (int jj = 0; jj < Field ::field_size; ++jj) {
                for (int i = 0; i < num_nodes; ++i) {
                    FaceBuild face = faces[i];
                    const EigMat<2, 1> x_f = face.getBarycenter();
                    const EigMat<2, 1> bdf = face.getBounds();
                    const EigMat<2, 2> rot = face.getRotationMatrix();
                    const EigMat<1, 1> s_f = face.getBarycenterF();
                    const EigMat<1, 1> bdf_proj = face.getBoundsF();
                    const real dist = (rot * (x_f - x_c))(1);
                    real normal_component_j;
                    if (dist > 0) {
                        normal_component_j = rot(1, jj);
                    } else {
                        normal_component_j = -rot(1, jj);
                    }
                    EigMat<cr, cl> m_adv_f = EigMat<cr, cl>::Zero();
                    EigMat<cr, fk> m_hyb_f = EigMat<cr, fk>::Zero();
                    for (int kk = 0; kk < FaceBuild ::dim_q; ++kk) {
                        const EigMat<2, 1> x_q_f = face.getQuadPoint(kk);
                        const EigMat<1, 1> s_q_f = face.getQuadPointF(kk);
                        const real w_q_f = face.getQuadWeight(kk);
                        const EigMat<cr, 1> d_phi_r_j = space.reco_basis.GetDerivativeVector(x_q_f, x_c, bdc, jj);
                        const EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                        const EigMat<fk, 1> psi_k = space.face_basis.GetEvaluationVector(s_q_f, s_f, bdf_proj);
                        m_adv_f += w_q_f * d_phi_r_j * phi_l.transpose();
                        m_hyb_f += w_q_f * d_phi_r_j * psi_k.transpose();
                    }
//                    EigMat<cr, fk> _lol = m_hyb_f * normal_component_j;
                    rhs.template block<cr, cl>(0, di * cl) -= m_adv_f * normal_component_j;
                    int fstart_i = field_size * cl + i * field_size * fk + di * fk;
                    rhs.template block<cr, fk>(0, fstart_i) += m_hyb_f * normal_component_j;
                }
            }
            return rhs.template block<cr - 1, elem_size>(1, 0);
        }

        EigMat<cr, elem_size> getCellProjectionOperator(intg di) const {
            EigMat<cr, elem_size> cell_projection_operator = EigMat<cr, elem_size>::Zero();
            EigMat<cl, cl> cell_eye_operator = EigMat<cl, cl>::Identity();
            cell_projection_operator.template block<cl, cl>(0, di * cl) += cell_eye_operator;
            return cell_projection_operator;
        }

        EigMat<cr * field_size, elem_size> getReconstructionOperator() const {
            EigMat<cr * field_size, elem_size> reconstruction_operator = EigMat<cr * field_size, elem_size>::Zero();
            const EigMat<cr - 1, cr - 1> lhs_inv = getReconstructionLHSInvert();
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            for (int i = 0; i < field_size; ++i) {
                EigMat<cr, elem_size> local_reconstruction_operator = EigMat<cr, elem_size>::Zero();
                EigMat<1, elem_size> line = EigMat<1, elem_size>::Zero();
                const EigMat<cr - 1, elem_size> rhs = getReconstructionRHS(i);
                local_reconstruction_operator.template block<cr - 1, elem_size>(1, 0) = lhs_inv * rhs;
                const EigMat<cr, elem_size> cell_projection_operator = getCellProjectionOperator(i);
                for (int j = 0; j < CellBuild ::dim_q; ++j) {
                    const EigMat<2, 1> x_q_c = cell.getQuadPoint(j);
                    const real w_q_c = cell.getQuadWeight(j);
                    const EigMat<cr, 1> phi_r = space.reco_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                    const real coef = (1.0 / real(CellBuild ::dim_q)) * w_q_c;
                    line += coef * phi_r.transpose() * (cell_projection_operator - local_reconstruction_operator);
                }
                reconstruction_operator.template block<cr - 1, elem_size>((i * cr) + 1, 0) += local_reconstruction_operator.template block<cr - 1, elem_size>(1, 0);
                reconstruction_operator.template block<1, elem_size>(i * cr, 0) += line;
            }
            return reconstruction_operator;
        }

        EigMat<elem_size, elem_size> getStabilizationOperator() const {
            EigMat<elem_size, elem_size> stabilization_operator_arg = EigMat<elem_size, elem_size>::Zero();
            EigMat<cl, cr> prj_r2l_t2t_rhs = EigMat<cl, cr>::Zero();
            EigMat<cl, cl> prj_cell_lhs = EigMat<cl, cl>::Zero();
            EigMat<cr * field_size, elem_size> recop = getReconstructionOperator();
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            for (int qc = 0; qc < CellBuild ::dim_q; ++qc) {
                const EigMat<2, 1> x_q_c = cell.getQuadPoint(qc);
                const real w_q_c = cell.getQuadWeight(qc);
                EigMat<cr, 1> phi_r = space.reco_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                prj_r2l_t2t_rhs += w_q_c * phi_l * phi_r.transpose();
                prj_cell_lhs += w_q_c * phi_l * phi_l.transpose();
            }
            EigMat<cl, cr> prj_r2l_t2t = prj_cell_lhs.llt().template solve(prj_r2l_t2t_rhs);
            for (int f = 0; f < num_nodes; ++f) {
                FaceBuild face = faces[f];
                const EigMat<2, 1> x_f = face.getBarycenter();
                const EigMat<2, 1> bdf = face.getBounds();
                const EigMat<2, 2> rot = face.getRotationMatrix();
                const EigMat<1, 1> s_f = face.getBarycenterF();
                const EigMat<1, 1> bdf_proj = face.getBoundsF();
                EigMat<fk, fk> proj_face_lhs = EigMat<fk, fk>::Zero();
                EigMat<fk, cr> prj_r2k_t2f_rhs = EigMat<fk, cr>::Zero();
                EigMat<fk, cl> prj_l2k_t2f_rhs = EigMat<fk, cl>::Zero();
                for (int qf = 0; qf < FaceBuild ::dim_q; ++qf) {
                    const EigMat<2, 1> x_q_f = face.getQuadPoint(qf);
                    const EigMat<1, 1> s_q_f = face.getQuadPointF(qf);
                    const real w_q_f = face.getQuadWeight(qf);
                    EigMat<fk, 1> psi_k = space.face_basis.GetEvaluationVector(s_q_f, s_f, bdf_proj);
                    EigMat<cr, 1> phi_r = space.reco_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    proj_face_lhs += w_q_f * psi_k * psi_k.transpose();
                    prj_r2k_t2f_rhs += w_q_f * psi_k * phi_r.transpose();
                    prj_l2k_t2f_rhs += w_q_f * psi_k * phi_l.transpose();
                }
                EigMat<fk, cr> prj_r2k_t2f = proj_face_lhs.llt().template solve(prj_r2k_t2f_rhs);
                EigMat<fk, cl> prj_l2k_t2f = proj_face_lhs.llt().template solve(prj_l2k_t2f_rhs);
                EigMat<fk, fk> face_eye = EigMat<fk, fk>::Identity();
                EigMat<cl, cl> cell_eye = EigMat<cl, cl>::Identity();
                for (int i = 0; i < field_size; ++i) {
                    EigMat<cl, elem_size> interop = EigMat<cl, elem_size>::Zero();
                    EigMat<fk, elem_size> stabilization_op33 = EigMat<fk, elem_size>::Zero();
                    interop.template block<cl, cl>(0, i * cl) += cell_eye;
                    EigMat<fk, elem_size> that_00 = prj_r2k_t2f * recop.template block<cr, elem_size>(i *cr, 0);
                    EigMat<cl, elem_size> that_0p = prj_r2l_t2t * recop.template block<cr, elem_size>(i *cr, 0);
                    EigMat<fk, elem_size> that_01 = prj_l2k_t2f * (- that_0p + interop);
                    EigMat<fk, elem_size> that = that_00 + that_01;
                    stabilization_op33 -= that;
                    const int fstart = field_size * cl + f * field_size * fk + i * fk;
                    stabilization_op33.template block<fk, fk>(0, fstart) += face_eye;
                    stabilization_operator_arg += (1.0 / bdf_proj.norm()) * stabilization_op33.transpose() * proj_face_lhs * stabilization_op33;
                }
            }
//            stabilization_operator_arg.print();
            return stabilization_operator_arg;
        }

        EigMat<ck, ck> getGradientLHSInvert() const {
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            EigMat<ck, ck> eye = EigMat<ck, ck>::Identity();
            EigMat<ck, ck> lhs = EigMat<ck, ck>::Zero();
            for (int i = 0; i < CellBuild::dim_q; ++i) {
                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
                const real w_q_c = cell.getQuadWeight(i);
                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                lhs += w_q_c * phik * phik.transpose();
            }
            EigMat<ck, ck> inv = lhs.llt().solve(eye);
            return inv;
        }

        EigMat<ck, elem_size> getGradientRHS(intg di, intg dj) const {
            EigMat<ck, elem_size> rhs = EigMat<ck, elem_size>::Zero();
            // CELL GEOMETRY
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
//            EigMat<ck, ck> lhs = EigMat<ck, ck>::Zero();
//            for (int i = 0; i < CellBuild::dim_q; ++i) {
//                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
//                const real w_q_c = cell.getQuadWeight(i);
//                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
//                lhs += w_q_c * phik * phik.transpose();
//            }
//            std::cout << lhs << std::endl;
            EigMat<ck, cl> m_adv_j = EigMat<ck, cl>::Zero();
            EigMat<ck, cl> m_adv_i = EigMat<ck, cl>::Zero();
            for (int i = 0; i < CellBuild::dim_q; ++i) {
                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
                const real w_q_c = cell.getQuadWeight(i);
                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                const EigMat<cl, 1> phil_i = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, di);
                const EigMat<cl, 1> phil_j = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, dj);
                m_adv_j += w_q_c * phik * phil_j.transpose();
                m_adv_i += w_q_c * phik * phil_i.transpose();
            }
//            std::cout << "m_adv_i" << std::endl;
//            std::cout << m_adv_i << std::endl;
//            std::cout << "m_adv_j" << std::endl;
//            std::cout << m_adv_j << std::endl;
            rhs.template block<ck, cl>(0, di * cl) += 1./2. * m_adv_j;
            rhs.template block<ck, cl>(0, dj * cl) += 1./2. * m_adv_i;
//            rhs.print();
            for (int i = 0; i < num_nodes; ++i) {
                FaceBuild face = faces[i];
                const EigMat<2, 1> x_f = face.getBarycenter();
                const EigMat<2, 1> bdf = face.getBounds();
                const EigMat<2, 2> rot = face.getRotationMatrix();
                const EigMat<1, 1> s_f = face.getBarycenterF();
                const EigMat<1, 1> bdf_proj = face.getBoundsF();
//                std::cout << "bdf_proj : " << bdf_proj << std::endl;
                const real dist = (rot * (x_f - x_c))(1);
                real normal_component_i;
                real normal_component_j;
                if (dist > 0) {
                    normal_component_i = rot(1, di);
                    normal_component_j = rot(1, dj);
                } else {
                    normal_component_i = -rot(1, di);
                    normal_component_j = -rot(1, dj);
                }
                EigMat<ck, cl> m_trace_f = EigMat<ck, cl>::Zero();
                EigMat<ck, fk> m_hyb_f = EigMat<ck, fk>::Zero();
                for (int j = 0; j < FaceBuild::dim_q; ++j) {
                    const EigMat<2, 1> x_q_f = face.getQuadPoint(j);
                    const EigMat<1, 1> s_q_f = face.getQuadPointF(j);
                    const real w_q_f = face.getQuadWeight(j);
                    const EigMat<ck, 1> phi_k = space.grad_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    const EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    const EigMat<fk, 1> psi_k = space.face_basis.GetEvaluationVector(s_q_f, s_f, bdf_proj);
                    m_trace_f += w_q_f * phi_k * phi_l.transpose();
                    m_hyb_f += w_q_f * phi_k * psi_k.transpose();
                }
                int fstart_i = field_size * cl + i * field_size * fk + di * fk;
                int fstart_j = field_size * cl + i * field_size * fk + dj * fk;
                rhs.template block<ck, cl>(0, di * cl) -= 1./2. * normal_component_j * m_trace_f;
                rhs.template block<ck, cl>(0, dj * cl) -= 1./2. * normal_component_i * m_trace_f;
                rhs.template block<ck, fk>(0, fstart_i) += 1./2. * normal_component_j * m_hyb_f;
                rhs.template block<ck, fk>(0, fstart_j) += 1./2. * normal_component_i * m_hyb_f;
            }
//            rhs.print();
//            real_matrix<ck, ck> eye = real_matrix<ck, ck>::Identity(ck, ck);
//            gr_lhs.ldlt().solve(gr_rhs);
//            EigMat<ck, elem_size> ret = lhs.llt().solve(rhs);
//            std::cout << ret << std::endl;
//            ret.print();
            return rhs;
        }

//        EigMat<ck, ck> getGradientLHSInvert() const {
//            const EigMat<2, 1> x_c = cell.getBarycenter();
//            const EigMat<2, 1> bdc = cell.getBounds();
//            EigMat<ck, ck> eye = EigMat<ck, ck>::Identity();
//            EigMat<ck, ck> lhs = EigMat<ck, ck>::Zero();
//            for (int i = 0; i < CellBuild::dim_q; ++i) {
//                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
//                const real w_q_c = cell.getQuadWeight(i);
//                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
//                lhs += w_q_c * phik * phik.transpose();
//            }
//            EigMat<ck, ck> inv = lhs.llt().solve(eye);
//            return inv;
//        }

        EigMat<ck, elem_size> getGradientFSRHS(intg di, intg dj) const {
            EigMat<ck, elem_size> rhs = EigMat<ck, elem_size>::Zero();
            // CELL GEOMETRY
            const EigMat<2, 1> x_c = cell.getBarycenter();
            const EigMat<2, 1> bdc = cell.getBounds();
            EigMat<ck, cl> m_adv_j = EigMat<ck, cl>::Zero();
//            EigMat<ck, cl> m_adv_i = EigMat<ck, cl>::Zero();
            for (int i = 0; i < CellBuild::dim_q; ++i) {
                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
                const real w_q_c = cell.getQuadWeight(i);
                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
//                const EigMat<cl, 1> phil_i = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, di);
                const EigMat<cl, 1> phil_j = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, dj);
                m_adv_j += w_q_c * phik * phil_j.transpose();
//                m_adv_i += w_q_c * phik * phil_i.transpose();
            }
            rhs.template block<ck, cl>(0, di * cl) += m_adv_j;
//            rhs.template block<ck, cl>(0, dj * cl) += 1./2. * m_adv_i;
            for (int i = 0; i < num_nodes; ++i) {
                FaceBuild face = faces[i];
                const EigMat<2, 1> x_f = face.getBarycenter();
                const EigMat<2, 1> bdf = face.getBounds();
                const EigMat<2, 2> rot = face.getRotationMatrix();
                const EigMat<1, 1> s_f = face.getBarycenterF();
                const EigMat<1, 1> bdf_proj = face.getBoundsF();
                const real dist = (rot * (x_f - x_c))(1);
//                real normal_component_i;
                real normal_component_j;
                if (dist > 0) {
//                    normal_component_i = rot(1, di);
                    normal_component_j = rot(1, dj);
                } else {
//                    normal_component_i = -rot(1, di);
                    normal_component_j = -rot(1, dj);
                }
                EigMat<ck, cl> m_trace_f = EigMat<ck, cl>::Zero();
                EigMat<ck, fk> m_hyb_f = EigMat<ck, fk>::Zero();
                for (int j = 0; j < FaceBuild::dim_q; ++j) {
                    const EigMat<2, 1> x_q_f = face.getQuadPoint(j);
                    const EigMat<1, 1> s_q_f = face.getQuadPointF(j);
                    const real w_q_f = face.getQuadWeight(j);
                    const EigMat<ck, 1> phi_k = space.grad_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    const EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_f, x_c, bdc);
                    const EigMat<fk, 1> psi_k = space.face_basis.GetEvaluationVector(s_q_f, s_f, bdf_proj);
                    m_trace_f += w_q_f * phi_k * phi_l.transpose();
                    m_hyb_f += w_q_f * phi_k * psi_k.transpose();
                }
                int fstart_i = field_size * cl + i * field_size * fk + di * fk;
//                int fstart_j = field_size * cl + i * field_size * fk + dj * fk;
                rhs.template block<ck, cl>(0, di * cl) -= normal_component_j * m_trace_f;
//                rhs.template block<ck, cl>(0, dj * cl) -= 1./2. * normal_component_i * m_trace_f;
                rhs.template block<ck, fk>(0, fstart_i) += normal_component_j * m_hyb_f;
//                rhs.template block<ck, fk>(0, fstart_j) += 1./2. * normal_component_i * m_hyb_f;
            }
            return rhs;
        }

        std::vector<EigMat<grad_size, elem_size>> getGradientOperators() const {
            std::vector<EigMat<grad_size, elem_size>> gradient_operators;
            gradient_operators.reserve(CellCmpt::dim_q);
            const EigMat<2, 1> x_c = cell_cmpt.getBarycenter();
            const EigMat<2, 1> bdc = cell_cmpt.getBounds();
            const EigMat<ck, ck> lhs_inv = getGradientLHSInvert();
            EigMat<grad_size, elem_size> gradient_operator = EigMat<grad_size, elem_size>::Zero();
            for (int i = 0; i < CellCmpt::dim_q; ++i) {
                const EigMat<2, 1> x_q_c = cell_cmpt.getQuadPoint(i);
                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
                for (int j = 0; j < Field::voigt_size; ++j) {
                    const intg row = std::get<0>(Field::voigt_data[j]);
                    const intg di = std::get<1>(Field::voigt_data[j]);
                    const intg dj = std::get<2>(Field::voigt_data[j]);
                    const real coef = std::get<3>(Field::voigt_data[j]);
                    EigMat<ck, elem_size> rhs;
                    if constexpr(Field::derivation_type == DerivationType::SYMMETRIC) {
                        rhs = getGradientRHS(di, dj);
                    }
                    else {
                        rhs = getGradientFSRHS(di, dj);
                    }
//                    const EigMat<ck, elem_size> rhs = getGradientRHS(di, dj);
                    const EigMat<ck, elem_size> grad_op = lhs_inv * rhs;
                    const EigMat<1, elem_size> line = phik.transpose() * grad_op;
//                    line.template block<1, elem_size>(0, 0).print();
//                    gradient_operator.template block<1, elem_size>(row, 0).print();
//                    gradient_operator.template block<1, elem_size>(row, 0) = coef * line;
                    gradient_operator.row(row) = coef * line;
                }
                gradient_operators.push_back(gradient_operator);
//                gradient_operator.print();
            }
            return gradient_operators;
        }

//        EigMat<ck, elem_size> getGradientComponent(intg di, intg dj) const {
//            EigMat<ck, elem_size> rhs = EigMat<ck, elem_size>::Zero();
//            // CELL GEOMETRY
//            const EigMat<2, 1> x_c = cell.getBarycenter();
//            const EigMat<2, 1> bdc = cell.getBounds();
//            EigMat<ck, ck> lhs = EigMat<ck, ck>::Zero();
//            for (int i = 0; i < CellBuild::dim_q; ++i) {
//                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
//                const real w_q_c = cell.getQuadWeight(i);
//                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
//                lhs += w_q_c * phik * phik.transpose();
//            }
//            std::cout << lhs << std::endl;
//            EigMat<ck, cl> m_adv_j = EigMat<ck, cl>::Zero();
//            EigMat<ck, cl> m_adv_i = EigMat<ck, cl>::Zero();
//            for (int i = 0; i < CellBuild::dim_q; ++i) {
//                const EigMat<2, 1> x_q_c = cell.getQuadPoint(i);
//                const real w_q_c = cell.getQuadWeight(i);
//                const EigMat<ck, 1> phik = space.grad_basis.GetEvaluationVector(x_q_c, x_c, bdc);
//                const EigMat<cl, 1> phil_i = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, di);
//                const EigMat<cl, 1> phil_j = space.cell_basis.GetDerivativeVector(x_q_c, x_c, bdc, dj);
//                m_adv_j += w_q_c * phik * phil_j.transpose();
//                m_adv_i += w_q_c * phik * phil_i.transpose();
//            }
//            std::cout << "m_adv_i" << std::endl;
//            std::cout << m_adv_i << std::endl;
//            std::cout << "m_adv_j" << std::endl;
//            std::cout << m_adv_j << std::endl;
//            rhs.template block<ck, cl>(0, di * cl) += 1./2. * m_adv_j;
//            rhs.template block<ck, cl>(0, dj * cl) += 1./2. * m_adv_i;
//            rhs.print();
//            for (int i = 0; i < num_nodes; ++i) {
//                FaceBuild face = faces[i];
//                const EigMat<2, 1> x_f = face.getBarycenter();
//                const EigMat<2, 1> bdf = face.getBounds();
//                const EigMat<2, 2> rot = face.getRotationMatrix();
//                const EigMat<1, 1> s_f = face.getBarycenterF();
//                const EigMat<1, 1> bdf_proj = face.getBoundsF();
////                std::cout << "bdf_proj : " << bdf_proj << std::endl;
//                const real dist = (rot * (x_f - x_c))(1);
//                real normal_component_i;
//                real normal_component_j;
//                if (dist > 0) {
//                    normal_component_i = rot(1, di);
//                    normal_component_j = rot(1, dj);
//                } else {
//                    normal_component_i = -rot(1, di);
//                    normal_component_j = -rot(1, dj);
//                }
//                EigMat<ck, cl> m_trace_f = EigMat<ck, cl>::Zero();
//                EigMat<ck, fk> m_hyb_f = EigMat<ck, fk>::Zero();
//                for (int j = 0; j < FaceBuild::dim_q; ++j) {
//                    const EigMat<2, 1> x_q_f = face.getQuadPoint(j);
//                    const EigMat<1, 1> s_q_f = face.getQuadPointF(j);
//                    const real w_q_f = face.getQuadWeight(j);
//                    const EigMat<ck, 1> phi_k = space.grad_basis.GetEvaluationVector(x_q_f, x_c, bdc);
//                    const EigMat<cl, 1> phi_l = space.cell_basis.GetEvaluationVector(x_q_f, x_c, bdc);
//                    const EigMat<fk, 1> psi_k = space.face_basis.GetEvaluationVector(s_q_f, s_f, bdf_proj);
//                    m_trace_f += w_q_f * phi_k * phi_l.transpose();
//                    m_hyb_f += w_q_f * phi_k * psi_k.transpose();
//                }
//                int fstart_i = field_size * cl + i * field_size * fk + di * fk;
//                int fstart_j = field_size * cl + i * field_size * fk + dj * fk;
//                rhs.template block<ck, cl>(0, di * cl) -= 1./2. * normal_component_j * m_trace_f;
//                rhs.template block<ck, cl>(0, dj * cl) -= 1./2. * normal_component_i * m_trace_f;
//                rhs.template block<ck, fk>(0, fstart_i) += 1./2. * normal_component_j * m_hyb_f;
//                rhs.template block<ck, fk>(0, fstart_j) += 1./2. * normal_component_i * m_hyb_f;
//            }
//            rhs.print();
////            real_matrix<ck, ck> eye = real_matrix<ck, ck>::Identity(ck, ck);
////            gr_lhs.ldlt().solve(gr_rhs);
//            EigMat<ck, elem_size> ret = lhs.llt().solve(rhs);
////            std::cout << ret << std::endl;
//            ret.print();
//            return ret;
//        }

    };
}

#endif //FETA_HYBRID_HXX
