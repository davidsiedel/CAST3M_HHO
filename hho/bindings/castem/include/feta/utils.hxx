//
// Created by dsiedel on 28/06/2021.
//

#ifndef FETA_UTILS_HXX
#define FETA_UTILS_HXX

namespace feta{
    static constexpr intg getBinomial(const intg n, const intg k) noexcept {
        return (k > n) ? 0 :                         // out of range
               (k == 0 || k == n) ? 1 :          // edge
               (k == 1 || k == n - 1) ? n :  // first
               (k + k < n)
               ?  // recursive:
               (getBinomial(n - 1, k - 1) * n) / k
               :  // path to face_polynomial_order=1   is faster
               (getBinomial(n - 1, k) * n) /
               (n - k);  // path to face_polynomial_order=n-1 is
        // faster
    }

    template <typename T, intg rows, intg cols>
    struct Array : std::array<T, rows * cols> {
        // default constructor
        constexpr Array() noexcept : std::array<T, rows * cols>() {}
        //        constexpr Array() noexcept = default;
        constexpr Array(Array &&) noexcept = default;
        constexpr Array(const Array &) noexcept = default;
        constexpr Array &operator=(Array &&) noexcept = default;
        constexpr Array &operator=(const Array &) noexcept = default;

        template <typename ValueType2>
        constexpr Array(const ValueType2 &v) {
            this->fill(v);
        }

        template <typename ValueType2>
        constexpr Array(const std::initializer_list<ValueType2> &values) {
            if ((values.size() != this->size()) && ((values.size() != 1))) {
                throw(std::runtime_error("invalid initializer list size"));
            }
            if (values.size() == 1) {
                this->fill(*(values.begin()));
            } else {
                auto p = this->begin();
                auto p2 = values.begin();
                for (; p != this->end(); ++p, ++p2) {
                    *p = *p2;
                }
            }
        }  // end of Array

        template <typename ValueType2>
        constexpr void fill(const ValueType2 &v) {
            for (auto p = this->begin(); p != this->end(); ++p) {
                *p = v;
            }
        }

        constexpr T &operator()(typename Array::size_type i,
                                typename Array::size_type j) {
            return std::array<T, rows * cols>::operator[](i * cols + j);
        }  // end operator()
        constexpr const T &operator()(typename Array::size_type i,
                                      typename Array::size_type j) const {
            return std::array<T, rows * cols>::operator[](i * cols + j);
        }  // end operator()

        static Array<T, rows, cols> get_transpose(const Array<T, rows, cols> &mat) {
            Array<T, rows, cols> temp;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    temp[j * cols + i] = mat[i * cols + j];
                }
            }
            return temp;
        }
    };
}

#endif //FETA_UTILS_HXX
