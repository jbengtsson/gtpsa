#ifndef _GTPSA_LIELIB_H_
#define _GTPSA_LIELIB_H_
#include <gtpsa/tpsa.hpp>
#include <gtpsa/ss_vect.h>
// #include <pybind11/stl.h>


class MNFType
{
private:
public:
  gtpsa::tpsa
    K,              // Lie generator in Floquet space.
    g;              /* Lie generator for canonical transformation to Floquet
		       space. */
  gtpsa::ss_vect<gtpsa::tpsa>
    M,              // Poincaré map.
    M_res,          // Residual map.
    A0, A0_inv,     // Linear transformation to fixed point.
    A1, A1_inv,     // Linear transformation to Floquet space.
    A_nl, A_nl_inv, // Nonlinear transformation to Floquet space.
    R;              // Floquet space rotation.

  MNFType(const std::shared_ptr<gtpsa::mad::desc> &desc, const int no):
    K(desc, no), g(desc, no), M(desc, no), M_res(desc, no), A0(desc, no),
    A0_inv(desc, no), A1(desc, no), A1_inv(desc, no), A_nl(desc, no),
    A_nl_inv(desc, no), R(desc, no)
  { }

};


void print_map(const std::string &str, const gtpsa::ss_vect<gtpsa::tpsa> &M);

void print_vec(const std::string &str, const std::vector<num_t> &v);

gtpsa::ss_vect<gtpsa::tpsa> param_to_ss_vect
(const int nm, const gtpsa::ss_vect<gtpsa::tpsa> &A,
 gtpsa::ss_vect<gtpsa::tpsa> &B);

gtpsa::ss_vect<gtpsa::tpsa> ss_vect_to_param
(const int nm, const gtpsa::ss_vect<gtpsa::tpsa> &A,
 gtpsa::ss_vect<gtpsa::tpsa> &B);
  
namespace gtpsa {
  /**
   * @brief Compute Dragt-Finn Lie generator from a (symplectic) map
   *
   */

  tpsa M_to_h(const ss_vect<tpsa> &M);
  tpsa M_to_h_DF(const ss_vect<tpsa> &t_map);
  ss_vect<tpsa> h_DF_to_M
  (const tpsa &h_DF, const ss_vect<tpsa> &x, const int k1, const int k2);

  void CtoR(const tpsa &a, tpsa &a_re, tpsa &a_im);
  tpsa RtoC(const tpsa &a_re, const tpsa &a_im);
  void GoFix(const ss_vect<tpsa> &map, ss_vect<tpsa> &A0);

} // namespace gtpsa

#endif //_GTPSA_LIELIB_H_
