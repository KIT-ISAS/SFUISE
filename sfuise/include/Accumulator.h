#pragma once

#include <array>
#include <chrono>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/CholmodSupport>

template <class T>
using SparseLLT = Eigen::CholmodSupernodalLLT<T>;

class SparseHashAccumulator
{
 public:

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorX;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
    typedef Eigen::Triplet<double> T;
    typedef Eigen::SparseMatrix<double> SparseMatrix;

    template <int ROWS, int COLS, typename Derived>
    inline void addH(int si, int sj, const Eigen::MatrixBase<Derived>& data)
    {
        EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, ROWS, COLS);
        KeyT id;
        id[0] = si;
        id[1] = sj;
        id[2] = ROWS;
        id[3] = COLS;
        auto it = hash_map.find(id);
        if (it == hash_map.end()) {
            hash_map.emplace(id, data);
        } else {
            it->second += data;
        }
    }

    template <int ROWS, typename Derived>
    inline void addB(int i, const Eigen::MatrixBase<Derived>& data)
    {
        b.template segment<ROWS>(i) += data;
    }

    inline void setup_solver()
    {
        std::vector<T> triplets;
        triplets.reserve(hash_map.size() * 36 + b.rows());
        for (const auto& kv : hash_map) {
            for (int i = 0; i < kv.second.rows(); i++) {
                for (int j = 0; j < kv.second.cols(); j++) {
                    triplets.emplace_back(kv.first[0] + i, kv.first[1] + j,
                    kv.second(i, j));
                }
            }
        }
        for (int i = 0; i < b.rows(); i++) {
            triplets.emplace_back(i, i, std::numeric_limits<double>::min());
        }
        smm = SparseMatrix(b.rows(), b.rows());
        smm.setFromTriplets(triplets.begin(), triplets.end());
    }

    inline VectorX Hdiagonal() const { return smm.diagonal(); }

    inline VectorX& getB() { return b; }

    inline VectorX solve(const VectorX* diagonal) const
    {
        auto t2 = std::chrono::high_resolution_clock::now();
        SparseMatrix sm = smm;
        if (diagonal) sm.diagonal() += *diagonal;
        VectorX res;
        SparseLLT<SparseMatrix> chol(sm);
        res = chol.solve(b);
        auto t3 = std::chrono::high_resolution_clock::now();
        auto elapsed2 = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
        return res;
    }

    inline void reset(int opt_size)
    {
        hash_map.clear();
        b.setZero(opt_size);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  private:

    using KeyT = std::array<int, 4>;

    struct KeyHash {
        inline size_t operator()(const KeyT& c) const
        {
            size_t seed = 0;
            for (int i = 0; i < 4; i++) {
                seed ^= c[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    std::unordered_map<KeyT, MatrixX, KeyHash> hash_map;
    VectorX b;
    SparseMatrix smm;
};
