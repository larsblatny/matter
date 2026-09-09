// Copyright (C) 2024 Lars Blatny. Released under GPL-3.0 license.

#ifndef SVD_HPP
#define SVD_HPP

#include <algorithm>
#include <limits>

#include "tools.hpp"

// Fast fixed-size SVD for the 2x2 / 3x3 deformation gradients used throughout the
// elastoplasticity kernel.
//
// Eigen::JacobiSVD is a general iterative solver; at these tiny fixed sizes it
// dominates plasticity(), which needs one full decomposition per particle per time
// step. This is a one-sided (Hestenes) Jacobi SVD specialized to DIMENSION:
// Givens rotations are applied to the columns of F from the right until the
// columns are mutually orthogonal. At that point
//
//     F V = W,     sigma_i = ||w_i||,     u_i = w_i / sigma_i,     F = U diag(sigma) V^T
//
// Why one-sided rather than the eigendecomposition of F^T F, which is the more
// obvious route at this size:
//
//   * Accuracy. Forming F^T F squares the condition number, so the small singular
//     values and, worse, their eigenvectors carry relative error O(eps*cond(F)^2).
//     One-sided Jacobi rotates the columns of F themselves by exact orthogonal
//     transforms and is relatively accurate in the sense of Demmel & Veselic:
//     error stays O(eps*cond(F)). This matters because the singular values feed
//     straight into the Hencky strain log(sigma).
//
//   * Orthogonality. V is a product of exact Givens rotations, so it is orthonormal
//     to machine precision however degenerate the singular values are. Closed-form
//     (Cardano) eigenvectors lose orthogonality precisely in the near-degenerate
//     case, which is the common one here: an F close to a rotation has all three
//     singular values equal.
//
// Conventions match Eigen::JacobiSVD with ComputeFullU | ComputeFullV, so this is a
// drop-in replacement at the call sites:
//   * F = U * diag(sigma) * V^T
//   * sigma >= 0, sorted in decreasing order
//   * U, V orthogonal, hence det(U) * det(V) = sign(det(F))
//
// Templated on the dimension rather than keyed off the DIMENSION macro, so that a
// build of either dimensionality can still instantiate and test both. FastSVD
// itself is the alias matching whatever TM/TV the rest of the code uses.
template <int D>
class FastSVDImpl {
public:
    using MatD = Eigen::Matrix<T, D, D>;
    using VecD = Eigen::Matrix<T, D, 1>;

    static_assert(D == 2 || D == 3, "FastSVD supports 2x2 and 3x3 only");

    explicit FastSVDImpl(const MatD& F) { compute(F); }

    const MatD& matrixU()        const { return U_; }
    const MatD& matrixV()        const { return V_; }
    const VecD& singularValues() const { return sigma_; }

private:
    MatD U_;
    MatD V_;
    VecD sigma_;

    void compute(const MatD& F) {
        // The sweeps below take dot products of columns, which squares magnitudes.
        // Prescaling by a power of two keeps those clear of overflow and underflow
        // for any |F| without perturbing a single digit: scaling by 2^-e and back
        // by 2^e is exact in binary floating point.
        // The shift is clamped so that neither the scaling 2^-exponent nor its undo
        // 2^exponent can overflow: frexp reports an exponent below -1023 for a
        // subnormal maxabs, and +1024 for a maxabs near the largest finite double,
        // and 2^1024 is not representable. Clamping still pulls both extremes far
        // enough in that the squared magnitudes stay in range, and powers of two
        // within this window are exact and normal, so the scaling stays lossless.
        static constexpr int max_shift = std::numeric_limits<T>::max_exponent - 24;

        const T maxabs = F.cwiseAbs().maxCoeff();
        int exponent = 0;
        if (maxabs > T(0) && std::isfinite(maxabs)) {
            std::frexp(maxabs, &exponent);
            exponent = std::max(-max_shift, std::min(max_shift, exponent));
        }

        MatD W = F * std::ldexp(T(1), -exponent); // columns, rotated in place
        V_.setIdentity();

        oneSidedJacobi(W);

        for (int i = 0; i < D; ++i)
            sigma_(i) = W.col(i).norm();

        sortDecreasing(W);
        buildU(W); // consumes the prescaled sigma_ for its rank threshold

        const T rescale = std::ldexp(T(1), exponent);
        for (int i = 0; i < D; ++i)
            sigma_(i) *= rescale;
    }

    // Rotate pairs of columns of W to mutual orthogonality, accumulating the
    // rotations into V_. Each pass over the column pairs reduces the off-diagonal
    // of the implicit Gram matrix quadratically; an F whose columns are already
    // orthogonal (F a rotation, or diagonal) costs one pass and no rotations.
    void oneSidedJacobi(MatD& W) {
        const T eps2 = std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon();
        const int max_sweeps = 16; // 2x2/3x3 converge in far fewer

        for (int sweep = 0; sweep < max_sweeps; ++sweep) {
            bool rotated = false;

            for (int p = 0; p < D - 1; ++p) {
                for (int q = p + 1; q < D; ++q) {
                    const T a = W.col(p).squaredNorm();
                    const T b = W.col(q).squaredNorm();
                    const T c = W.col(p).dot(W.col(q));

                    // Skip if the pair is already orthogonal to working precision.
                    // Negated so a NaN also takes this branch, and so that tau
                    // below is bounded, keeping tau*tau clear of overflow.
                    if (!(c * c > eps2 * a * b))
                        continue;

                    // Rotation diagonalizing the 2x2 Gram matrix [[a, c], [c, b]];
                    // the root of t^2 + 2*tau*t - 1 = 0 of smaller magnitude is
                    // taken, which keeps the rotation angle within +/- pi/4.
                    const T tau = (b - a) / (T(2) * c);
                    const T w   = std::sqrt(tau * tau + T(1));
                    const T t   = (tau >= T(0)) ? T(1) / (tau + w) : T(1) / (tau - w);
                    const T cs  = T(1) / std::sqrt(t * t + T(1));
                    const T sn  = t * cs;

                    for (int r = 0; r < D; ++r) {
                        const T wp = W(r, p);
                        const T wq = W(r, q);
                        W(r, p) = cs * wp - sn * wq;
                        W(r, q) = sn * wp + cs * wq;

                        const T vp = V_(r, p);
                        const T vq = V_(r, q);
                        V_(r, p) = cs * vp - sn * vq;
                        V_(r, q) = sn * vp + cs * vq;
                    }
                    rotated = true;
                }
            }

            if (!rotated)
                break;
        }
    }

    // Sort the singular values into decreasing order, permuting the columns of V_
    // and W with them. Each swap flips det(V_), which is exactly what preserves
    // det(U)*det(V) = sign(det(F)) once U is built from the permuted W.
    void sortDecreasing(MatD& W) {
        auto swapIfNeeded = [&](int i, int j) {
            if (sigma_(i) < sigma_(j)) {
                std::swap(sigma_(i), sigma_(j));
                V_.col(i).swap(V_.col(j));
                W.col(i).swap(W.col(j));
            }
        };
        // Sorting networks for the only two sizes this supports.
        if constexpr (D == 3) {
            swapIfNeeded(0, 1);
            swapIfNeeded(1, 2);
            swapIfNeeded(0, 1);
        } else {
            swapIfNeeded(0, 1);
        }
    }

    // Normalize the columns of W into U_. They are already orthogonal to working
    // precision, so the Gram-Schmidt below is nearly a no-op in the ordinary case;
    // it costs a handful of flops and bounds ||U^T U - I|| when F is ill
    // conditioned. A column whose norm has collapsed carries no direction
    // information (F is rank deficient there, or sigma_i has sunk to the noise
    // floor of the product F*V) and is replaced by an orthonormal complement.
    //
    // The projection is subtracted twice. One pass leaves O(eps * cond(F)) loss of
    // orthogonality, because for a column nearly dependent on the earlier ones the
    // subtraction cancels and the small residual keeps a large relative error; a
    // second pass restores O(eps) ("twice is enough" -- Kahan and Parlett). The
    // difference is invisible in double but reaches 4e-5 in float.
    void buildU(const MatD& W) {
        const T tol = T(16) * std::numeric_limits<T>::epsilon() * sigma_(0);

        for (int i = 0; i < D; ++i) {
            VecD u = W.col(i);
            for (int pass = 0; pass < 2; ++pass)
                for (int j = 0; j < i; ++j)
                    u -= U_.col(j).dot(u) * U_.col(j);

            const T n = u.norm();
            U_.col(i) = (n > tol) ? VecD(u / n) : complement(i);
        }
    }

    // A unit vector orthogonal to the first i columns of U_.
    VecD complement(int i) const {
        if constexpr (D == 3) {
            if (i == 0)
                return VecD(1, 0, 0);
            if (i == 1) {
                // Orthogonalize against whichever axis the column we have leans on
                // least; that axis is at most 1/sqrt(3) aligned with it, so the
                // residual norm is at least 0.8 and this is never near degenerate.
                // Projecting the axis rather than crossing with it keeps the result
                // canonical: an axis-aligned U_.col(0) yields the next axis, so a
                // zero F comes back as the identity, not a permutation of it.
                const VecD u0 = U_.col(0);
                const T ax = std::abs(u0(0)), ay = std::abs(u0(1)), az = std::abs(u0(2));
                const VecD axis = (ax <= ay && ax <= az) ? VecD(1, 0, 0)
                                                         : ((ay <= az) ? VecD(0, 1, 0)
                                                                       : VecD(0, 0, 1));
                return VecD(axis - u0.dot(axis) * u0).normalized();
            }
            return U_.col(0).cross(U_.col(1));
        } else {
            if (i == 0)
                return VecD(1, 0);
            return VecD(-U_(1, 0), U_(0, 0));
        }
    }
};

// The decomposition at the dimensionality the rest of the code is compiled for.
using FastSVD = FastSVDImpl<DIMENSION>;

#endif // SVD_HPP
