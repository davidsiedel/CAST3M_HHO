//
// Created by dsiedel on 28/04/2021.
//

#ifndef FETA_FETA_HXX
#define FETA_FETA_HXX

//#include <Eigen/Dense>
#include "lib/eigen3/Eigen/Dense"
#include <iostream>
#include "config/config.hxx"

namespace feta {

//    using namespace Eigen;

    template<intg r, intg c>
    using EigMap = Eigen::Map<const Eigen::Matrix<real, r, c, Eigen::RowMajor>>;
    template<intg r, intg c>
    using EigMapC = Eigen::Map<const Eigen::Matrix<real, r, c, Eigen::ColMajor>>;
    template<intg r, intg c>
    using EigMapIntR = Eigen::Map<const Eigen::Matrix<intg, r, c, Eigen::RowMajor>>;
    template<intg r, intg c>
    using EigMapIntC = Eigen::Map<const Eigen::Matrix<intg, r, c, Eigen::ColMajor>>;

//    template<intg r, intg c>
//    using EigMat = Eigen::Matrix<real, r, c>;
    template<intg r, intg c>
    using EigMat2 = Eigen::Matrix<intg, r, c>;

    template<intg r, intg c>
    struct EigMat : public Eigen::Matrix<real, r, c>{
        using MAT = Eigen::Matrix<real, r, c>;
        using MAT::MAT;

        void print() const {
//            Eigen::IOFormat PrintFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
            Eigen::IOFormat PrintFmt(4, 0, ", ", ",\n", "[", "]", "[", "]");
            std::cout << std::fixed << this->format(PrintFmt) << std::endl;
        }

    };

    template<typename T>
    struct GenericMatrix : public T{
        using T::T;
        void print(){
//            Eigen::IOFormat PrintFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
            Eigen::IOFormat PrintFmt(4, 0, ", ", ",\n", "[", "]", "[", "]");
            std::cout << std::fixed << this->format(PrintFmt) << std::endl;
        }
    };

    template<intg r, intg c>
    using Mat = GenericMatrix<Eigen::Matrix<real, r, c>>;

    template<intg c>
    using Vec = GenericMatrix<Eigen::Matrix<real, 1, c>>;

    enum struct QuadratureType {
        Gauss
    };

    enum struct BasisType {
        Monomial
    };

    static constexpr real SQRT2 = 1.4142135623730950488016887242096980785696718753;

    template<typename matrix_t>
    void print(matrix_t mat){
        Eigen::IOFormat PrintFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
        std::cout << std::fixed << mat.format(PrintFmt) << std::endl;
    }

}

#endif //FETA_FETA_HXX
