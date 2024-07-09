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
    K,     // Lie generator in Floquet space.
    g;     // Lie generator for canonical transformation to Floquet space.
  gtpsa::ss_vect<gtpsa::tpsa>
    M,     // Poincaré map.
    A_0,   // Linear transformation to fixed point.
    A_1,   // Linear transformation to Floquet space.
    A_nl,  // Nonlinear transformation to Floquet space.
    R;     // Map in Floquet space.

  MNFType(const std::shared_ptr<gtpsa::mad::desc> &desc, const int no):
    K(desc, no), g(desc, no), M(desc, no), A_0(desc, no), A_1(desc, no),
    A_nl(desc, no), R(desc, no)
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
  tpsa M_to_h_DF(const ss_vect<tpsa> &M);
  ss_vect<tpsa> h_DF_to_M
  (const tpsa &h_DF, const ss_vect<tpsa> &x, const int k1, const int k2);

  void CtoR(const tpsa &a, tpsa &a_re, tpsa &a_im);
  tpsa RtoC(const tpsa &a_re, const tpsa &a_im);
  ss_vect<tpsa> GoFix(const ss_vect<tpsa> &M);
  MNFType Map_Norm(const ss_vect<tpsa> &M);

} // namespace gtpsa

#endif //_GTPSA_LIELIB_H_
