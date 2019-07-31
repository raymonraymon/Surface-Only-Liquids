#pragma once
/** @file LaplaceSpherical.hpp
 * @brief Implements the Laplace kernel with spherical expansions.
 *
 * K(t,s) = 1 / |s-t|        // Laplace potential
 * K(t,s) = (s-t) / |s-t|^3  // Laplace force
 */

#include <complex>
#include <cmath>
#include <cassert>

#include "fmmtl/Expansion.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"

#include "kernel/Util/SphericalMultipole3D.hpp"

#include "Laplace.kern"
#include "../../Timer.h"

class LaplaceSpherical
    : public fmmtl::Expansion<LaplaceKernel, LaplaceSpherical> {
 public:
  typedef double real_type;
  typedef std::complex<real_type> complex_type;

  //! Point type
  typedef Vec<3,real_type> point_type;

  //! Multipole expansion type
  typedef std::vector<complex_type> multipole_type;
  //! Local expansion type
  typedef std::vector<complex_type> local_type;

  //! Transform operators
  typedef SphericalMultipole3D<point_type,multipole_type,local_type> SphOp;

  //! Expansion order
  int P;

  //! Constructor
  LaplaceSpherical(int _P = 5)
      : P(_P) {
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M = std::vector<complex_type>(P*(P+1)/2);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, const point_type&, unsigned) const {
    L = std::vector<complex_type>(P*(P+1)/2);
  }

  /** Kernel S2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void S2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
//      CSim::TimerMan::timer("FMM_multiply/execute/upward/SLP_S2M").start();
    SphOp::S2M(P, center-source, charge, M);
//      CSim::TimerMan::timer("FMM_multiply/execute/upward/SLP_S2M").stop();
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all points within its box
   */
  void M2M(const multipole_type& Msource,
           multipole_type& Mtarget,
           const point_type& translation) const {
//      CSim::TimerMan::timer("FMM_multiply/execute/upward/SLP_M2M").start();
    SphOp::M2M(P, Msource, Mtarget, translation);
//      CSim::TimerMan::timer("FMM_multiply/execute/upward/SLP_M2M").stop();
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] Msource The multpole expansion source
   * @param[in,out] Ltarget The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre Msource includes the influence of all points within its box
   */
  void M2L(const multipole_type& Msource,
           local_type& Ltarget,
           const point_type& translation) const {
//      CSim::TimerMan::timer("FMM_multiply/execute/far/SLP_M2L").start();
    SphOp::M2L(P, Msource, Ltarget, translation);
//      CSim::TimerMan::timer("FMM_multiply/execute/far/SLP_M2L").stop();
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& Lsource,
           local_type& Ltarget,
           const point_type& translation) const {
//      CSim::TimerMan::timer("FMM_multiply/execute/downward/SLP_L2L").start();
    SphOp::L2L(P, Lsource, Ltarget, translation);
//      CSim::TimerMan::timer("FMM_multiply/execute/downward/SLP_L2L").stop();
  }

  /** Kernel L2T operation
   * r += Op(L, t) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] target The target of this L2T operation
   * @param[in] result The result to accumulate into
   * @pre L includes the influence of all sources outside its box
   */
  void L2T(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
//      CSim::TimerMan::timer("FMM_multiply/execute/downward/SLP_L2T").start();
    using std::real;
    using std::imag;

    real_type rho, theta, phi;
    SphOp::cart2sph(rho, theta, phi, target - center);
      if (theta == 0) theta += 1e-100;

    complex_type Z[P*(P+1)/2], dZ[P*(P+1)/2];
    SphOp::evalZ(rho, theta, phi, P, Z, dZ);

    point_type sph = point_type();
    int nm = 0;
    for (int n = 0; n != P; ++n) {
      const real_type LZ = real(L[nm])*real(Z[nm]) - imag(L[nm])*imag(Z[nm]);
      result[0] += LZ;
      sph[0]    += LZ / rho * n;
      sph[1]    += real(L[nm])*real(dZ[nm]) - imag(L[nm])*imag(dZ[nm]);

      ++nm;
      for (int m = 1; m <= n; ++m, ++nm) {
        const complex_type LZ = L[nm] * Z[nm];
        result[0] += 2 * real(LZ);
        sph[0]    += 2 * real(LZ) / rho * n;
        sph[1]    += 2 * (real(L[nm])*real(dZ[nm])-imag(L[nm])*imag(dZ[nm]));
        sph[2]    += 2 *-imag(LZ) * m;
      }
    }

    const point_type cart = SphOp::sph2cart(rho, theta, phi, sph);
    result[1] += cart[0];
    result[2] += cart[1];
    result[3] += cart[2];
//      CSim::TimerMan::timer("FMM_multiply/execute/downward/SLP_L2T").stop();
  }
};
