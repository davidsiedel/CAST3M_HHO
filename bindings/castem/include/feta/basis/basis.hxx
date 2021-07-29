//
// Created by dsiedel on 28/06/2021.
//

#ifndef FETA_BASIS_HXX
#define FETA_BASIS_HXX

#include "feta/feta.hxx"
#include "feta/utils.hxx"

namespace feta::basis {

    template<intg d, intg ord>
    struct Monomial {
        Monomial() : exponents(getExps()) {}

        static constexpr intg dim = getBinomial(d + ord, d);
        static constexpr intg order = ord;


        EigMat<dim, 1>
        GetEvaluationVector(const EigMat<d, 1> &point, const EigMat<d, 1> &centroid, const EigMat<d, 1> &bounds) const {
            EigMat<dim, 1> v;
            for (int i = 0; i < dim; ++i) {
                real value = 1.0;
                for (int j = 0; j < d; ++j) {
                    value *= std::pow(2.0 * (point(j) - centroid(j)) / bounds(j), exponents(i, j));
                }
                v(i) = value;
            }
            return v;
        }


        EigMat<dim, 1>
        GetDerivativeVector(const EigMat<d, 1> &point, const EigMat<d, 1> &centroid, const EigMat<d, 1> &bounds,
                            intg dx) const {
            EigMat<dim, 1> v;
            for (int i = 0; i < dim; ++i) {
                real value = 1.0;
                for (int j = 0; j < d; ++j) {
                    if (j != dx) {
                        value *= std::pow(2.0 * (point(j) - centroid(j)) / bounds(j), exponents(i, j));
                    } else {
                        if (exponents(i, j) > 0) {
                            real coef = 2.0 * exponents(i, j) / bounds(j);
                            value *= coef * std::pow(2.0 * (point(j) - centroid(j)) / bounds(j), exponents(i, j) - 1);
                        } else {
                            value *= 0.0;
                        }
                    }
                }
                v(i) = value;
            }
            return v;
        }

        EigMat2<dim, d> exponents;

    private:
        EigMat2<dim, d> getExps() const {
            EigMat2<dim, d> a;
            if constexpr (d == 1) {
//            if (d == 1) {
                int row_count = 0;
                for (int i = 0; i < ord + 1; ++i) {
                    a(row_count, 0) = intg(i);
                    row_count += 1;
                }
            } else if constexpr (d == 2) {
//            } else if (d == 2) {
                int row_count = 0;
                for (int i = 0; i < ord + 1; ++i) {
                    for (int j = 0; j < i + 1; ++j) {
                        int m = i - j;
                        a(row_count, 0) = intg(m);
                        a(row_count, 1) = intg(j);
                        row_count += 1;
                    }
                }
            } else if constexpr (d == 3) {
//            } else if (d == 3) {
                int row_count = 0;
                for (int i = 0; i < ord + 1; ++i) {
                    for (int j = 0; j < i + 1; ++j) {
                        for (int k = 0; k < i + 1; ++k) {
                            if (j + k < i + 1) {
                                int m = i - (j + k);
                                a(row_count, 0) = intg(m);
                                a(row_count, 1) = intg(k);
                                a(row_count, 2) = intg(j);
                                row_count += 1;
                            }
                        }
                    }
                }
            }
            return a;
        }
    };

}

#endif //FETA_BASIS_HXX
