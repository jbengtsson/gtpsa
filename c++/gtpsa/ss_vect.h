#ifndef _GTPSA_SS_VECT_H_
#define _GTPSA_SS_VECT_H_ 1

#include <gtpsa/desc.hpp>
#include <gtpsa/tpsa.hpp>
//#include <gtpsa/intern/gtpsa_container.hpp>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <armadillo>

// #include <gtpsa/lielib.hpp>

namespace gtpsa {

  const int ss_vect_n_dim = 6 + 1;


  /**
   * state space, phase space vector
   *
   * assuming that all are of the same length ...
   *
   * @todo:
   *        review space state length
   *        review description the state space should have
   */

  void  report_vector_dimension_mismatch(size_t nv,  const size_t n);

  inline void check_vector_dimension
  (const std::shared_ptr<mad::desc> desc, const size_t n)
  {
    size_t nv = size_t(desc->getNv());
    if (n > nv){
      report_vector_dimension_mismatch(nv, n);
    }
  }


  template<typename T>
  class ss_vect {
  private:
    std::vector<T> state_space;
    std::string m_name = "";

  public:
    inline ss_vect
    (const std::shared_ptr<mad::desc> desc, const ord_t mo,
     const size_t n = ss_vect_n_dim)
    {
#if 1
      check_vector_dimension(desc, n);
#endif
      this->state_space.reserve(n);
      for (size_t i = 0; i < n; ++i) {
	T tmp(desc, mo);
	this->state_space.push_back(tmp);
      }
    }

    inline ss_vect(const T &t, const size_t n = ss_vect_n_dim) {
#if 1
      check_vector_dimension(t.getDescription(), n);
#endif
      this->state_space.reserve(n);
      for (size_t i = 0; i < n; ++i) {
	this->state_space.push_back(t.newFromThis());
      }

    }


    /*tps
      inline ss_vect(const T& t, std::vector<double> vec){
      this->state_space.reserve(vec.size());
      for(size_t i=0; i<vec.size(); ++i){
      this->state_space.push_back(vec[i]);
      }
      }
    */

    inline ss_vect
    (const T &x, const T &px, const T &y, const T &py, const T &delta,
     const T &ct)
      : state_space{x, px, y, py, delta, ct}
    { this->checkSize(this->state_space); }

    /*
      this->state_space.reserve(6);
      this->state_space.push_back(x);
      this->state_space.push_back(px);
      this->state_space.push_back(y);
      this->state_space.push_back(py);
      this->state_space.push_back(delta);
      this->state_space.push_back(ct);
      }
    */

    inline ss_vect(const std::vector<T> &vec)
      : state_space(vec)
    { this->checkSize(this->state_space); }

    void setName(std::string name) {
      this->m_name = name;
    }

    std::string name(void) {
      return this->m_name;
    }

    inline ss_vect<T> &operator+=(const ss_vect<T> &o) {
      this->inplace_op(o, [](T &t, const T &o) { t += o; });
      return *this;
    }

    inline ss_vect<T> &operator-=(const ss_vect<T> &o) {
      this->inplace_op(o, [](T &t, const T &o) { t -= o; });
      return *this;
    }

    /*
    inline ss_vect<T>&  operator*=
    (const ss_vect<T>& o){ this->inplace_op(o, [] (T& t, const T&     o)
    { t *= o;}); return *this; }
    inline ss_vect<T>&  operator/=
    (const ss_vect<T>& o){ this->inplace_op(o, [] (T& t, const T&     o)
    { t /= o;}); return *this; }
    */


    inline ss_vect<T> &operator+=(const double o) {
      this->inplace_op(o, [](T &t, const double o) { t += o; });
      return *this;
    }

    inline ss_vect<T> &operator-=(const double o) {
      this->inplace_op(o, [](T &t, const double o) { t -= o; });
      return *this;
    }

    inline ss_vect<T> &operator*=(const double o) {
      this->inplace_op(o, [](T &t, const double o) { t *= o; });
      return *this;
    }

    inline ss_vect<T> &operator/=(const double o) {
      this->inplace_op(o, [](T &t, const double o) { t /= o; });
      return *this;
    }

    template<typename U>
    inline ss_vect<T> &operator+=(const ss_vect<U> &o) {
      this->inplace_op<U>(o, [](T &t, const U &o) { t += o; });
      return *this;
    }

    template<typename U>
    inline ss_vect<T> &operator-=(const ss_vect<U> &o) {
      this->inplace_op<U>(o, [](T &t, const U &o) { t -= o; });
      return *this;
    }

    /*
    template<typename U>
    inline ss_vect<T>& operator*=
    (const ss_vect<double>& o){ this->inplace_op(o, [] (T& t, const T&     o)
    { t *= o;}); return *this; }
    template<typename U>
    inline ss_vect<T>& operator/=
    (const ss_vect<double>& o){ this->inplace_op(o, [] (T& t, const T&     o)
    { t /= o;}); return *this; }
    */

    inline void set_zero(void);
    /*{
      std::for_each(this->state_space.begin(), this->state_space.end(),
		    [](gtpsa::tpsa& v){v = 0.0;});
    }*/

    inline void set_identity(void) {
      throw std::runtime_error("set_identity only implemented for tpsa");
      //
    }

    inline auto size(void) const { return this->state_space.size(); }

#ifndef GTSPA_ONLY_OPTIMISED_OPS

    inline ss_vect<T> operator+(const ss_vect<T> &o) const {
      auto r = this->clone();
      r += o;
      return r;
    }

    inline ss_vect<T> operator-(const ss_vect<T> &o) const {
      auto r = this->clone();
      r -= o;
      return r;
    }

    /*
    inline ss_vect<T> operator * (const ss_vect<T>& o) const
    { auto r = this->clone(); r *= o; return r; }
    inline ss_vect<T> operator / (const ss_vect<T>& o) const
    { auto r = this->clone(); r /= o; return r; }
    */

    inline ss_vect<T> operator+(const double o) const {
      auto r = this->clone();
      r += o;
      return r;
    }

    inline ss_vect<T> operator-(const double o) const {
      auto r = this->clone();
      r -= o;
      return r;
    }

    inline ss_vect<T> operator*(const double o) const {
      auto r = this->clone();
      r *= o;
      return r;
    }

    inline ss_vect<T> operator/(const double o) const {
      auto r = this->clone();
      r /= o;
      return r;
    }

    /*
    ss_vect(const ss_vect& o) {
      this->state_space.reserve(o.state_space.size());
      std::copy
	(o.state_space.begin(), o.state_space.end(),
	 std::back_inserter(this->state_space));
    }
    */
    ss_vect(const ss_vect &o) = default;

#else
    ss_vect(const ss_vect& o) = delete;
#endif /* GTSPA_ONLY_OPTIMISED_OPS */

    inline ss_vect(const ss_vect<T> &&o)
      noexcept: state_space(std::move(o.state_space)) {};

    inline ss_vect<T> &operator=(const ss_vect<T> &&o) noexcept {
      if (this == &o) {
	return *this;
      }
      this->state_space = std::move(o.state_space);
      // o->state_space = nullptr;
      return *this;
    }

    /**
     * @brief allocates empty objects with the same property as this object
     *
     * similar to clone, but it does not copy the content of this object
     *
     *  allocateLikeMe ?
     */
    inline ss_vect<T> allocateLikeMe(void) const {
      const std::vector<T> &vec = this->state_space;
      this->checkSize(vec);

      // back inserted as the single argument is not default constructable
      std::vector<T> nvec;
      nvec.reserve(this->state_space.size());
      std::transform(vec.begin(), vec.end(), std::back_inserter(nvec),
		     [](const T &elem) -> T { return same_as_instance(elem); });
      auto ret =  ss_vect(nvec);
      // a paranoic extra
      ret.set_zero();
      return ret;
    }
#if 0
    /* seems the other name won the race */
    inline ss_vect<T> newFromThis(void) const {
      throw std::runtime_error("Implement me but decide on method name first!");
    }
#endif
    /**
     * @todo implement a more flexible clone
     */
    inline ss_vect<T> clone(void) const {
      const std::vector<T> &vec = this->state_space;
      this->checkSize(vec);

      // back inserted as the single argument is not default constructable
      std::vector<T> nvec;
      nvec.reserve(this->state_space.size());
      std::transform(vec.begin(), vec.end(), std::back_inserter(nvec),
		     [](const T &elem) -> T { return elem.clone(); });

      return ss_vect(nvec);
    }

    inline void _copyInPlace(const ss_vect<T> &o) {
      std::transform
	(o.state_space.begin(), o.state_space.end(), this->state_space.begin(),
	 [](const T &elem) -> T { return elem.clone(); });
    }

    inline T &operator[](const int i) { return this->state_space[i]; }

    inline const T &operator[](const int i) const
    { return this->state_space[i]; }

    inline T &at(const int i) { return this->state_space.at(i); }

    inline const T &at(const int i) const { return this->state_space.at(i); }

    inline const std::vector<T> getConstVector(void) const
    {return this->state_space;}
    inline auto getVector(void) const { return this->state_space; }
    inline auto getVector(void)       { return this->state_space; }

    inline void cst(std::vector<double> &r) const {
      for (size_t i = 0; i < this->size(); ++i)
	{ r[i] = this->state_space[i].cst(); }
      // std::transform(begin(), this->state_space.end(), r.state_space.begin(),
      // 		     [](T& elem) { return elem.cst();});
    }

    inline std::vector<double> cst(void) const {
      std::vector<double> r(this->state_space.size());
      this->cst(r);
      return r;
    }


    void show(std::ostream &strm, int level = 1, bool with_endl = true) const;

    std::string pstr(void);

    std::string repr(void);


    inline arma::mat toMatrix(void) {
      throw std::runtime_error("only implemented for tps(a)");
    }

    inline arma::mat jacobian(void) const {
      throw std::runtime_error("jacobian only implemented for tps(a)");
    }

    inline void setJacobian(arma::mat &jac) {
      throw std::runtime_error("setJacobian only implemented for tps(a)");
    }

    inline arma::cube hessian() const {
      throw std::runtime_error("hessian only implemented for tps(a)");
    }

    inline void setHessian(arma::cube& jac) {
      throw std::runtime_error("setHessian only implemented for tps(a)");
    }

    inline void setConstantPartWithoutCheck(const std::vector<double> &vec) {
      throw std::runtime_error("Needs to be implemented as specialisation");
    }

    inline void setConstant(const std::vector<double> &vec) {
      if (vec.size() != this->state_space.size()) {
	std::stringstream strm;
	strm << "Vector size " << vec.size()
	     << " does not match size of space state "
	     << this->state_space.size();
	throw std::runtime_error(strm.str());
      }
      this->setConstantPartWithoutCheck(vec);
    }

    inline void rminv(const ss_vect<T>& a) {
      throw std::runtime_error("rminv currently only implemented for tpsa");
    }

    inline void rpminv(const ss_vect<T>& a, std::vector<idx_t>& select) {
      throw std::runtime_error("rpminv currently only implemented for tpsa");
    }

    inline void rcompose(const ss_vect<T>& a, const ss_vect<T>& b) {
      throw std::runtime_error("rcompose currently only implemented for tpsa");
    }

    void rvec2fld(const T& a) {
      throw std::runtime_error("rvec2fld currently only implemented for tpsa");
    }

    void fld2vec(T* a) const {
      throw std::runtime_error("fld2vec currently only implemented for tpsa");
    }

    void fgrad(T* a, T*b) const {
      throw std::runtime_error("fgrad currently only implemented for tpsa");
    }

    void rliebra(const ss_vect<T>& a, const ss_vect<T>& b) {
      throw std::runtime_error("rliebra currently only implemented for tpsa");
    }

    void rexppb(const ss_vect<T>& a, const ss_vect<T>& b) {
      throw std::runtime_error("rexpbp currently only implemented for tpsa");
    }

    void rlogpb(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b) {
      throw std::runtime_error("rlogpb currently only implemented for tpsa");
    }

    void rgetOrder(const ss_vect<T>& o, const int order) {
      throw std::runtime_error("currently only implemented for tpsa");
    }

    void rderiv(const ss_vect<T>& o, const int order) {
      throw std::runtime_error("currently only implemented for tpsa");
    }


    void get_mns(const int no1, const int no2, ss_vect<tpsa> &y) const;

    void M_to_h(tpsa &h) const;

    void M_to_h_DF(tpsa &h) const;

    void GoFix(ss_vect<tpsa> &A_0) const;

    void Map_Norm
    (ss_vect<tpsa> &A_0, ss_vect<tpsa> &A_1, ss_vect<tpsa> &R, tpsa &g, tpsa &K)
      const;


    /**
     * @btodo:
     *    default to zero?
     * */
    inline int getMaximumOrder(void) const  {
      throw std::runtime_error
	("maximum order currently only implemented for tpsa");
    }

    inline double computeNorm(void) const {
      throw std::runtime_error("norm to be provided for xxx");
    }

  private:
    inline void checkSize(const std::vector<T> &vec) const {
      if (vec.size() != this->size()) {
	std::ostringstream strm;
	strm << "Currently ss_vect implemented for length " << this->size()
	     <<  " got length  " << vec.size();
	throw std::runtime_error(strm.str());
      }
    }

    // I am sure there it is a way to use std algorithms ...
    inline void inplace_op(const ss_vect<T> &o, void f(T &, const T &)) {
      for (size_t i = 0; i < this->state_space.size(); ++i) {
	f(this->state_space[i], o.state_space[i]);
      }
    }

    template<typename U>
    inline void inplace_op(const ss_vect<U> &o, void f(T &, const U &)) {
      for (size_t i = 0; i < this->state_space.size(); ++i) {
	f(this->state_space[i], o[i]);
      }
    }

    inline void inplace_op(const double o, void f(T &, const double)) {
      for (size_t i = 0; i < this->state_space.size(); ++i) {
	f(this->state_space[i], o);
      }
    }
  };

  template<typename T>
  inline std::ostream &operator<<(std::ostream &strm, const ss_vect<T> &s) {
    s.show(strm, 1);
    return strm;
  }


  /*
   * @brief special installation
   */
  template<>
  inline ss_vect<double>::ss_vect
  (const std::shared_ptr<mad::desc> d, const ord_t m, const size_t n)
    : state_space(n)
  {}

  template<>
  inline ss_vect<double>::ss_vect(const double& d, const size_t n)
    : state_space(n)
  {}

  template<>
  inline ss_vect<std::complex<double>>::ss_vect
  (const std::shared_ptr<mad::desc> d, const ord_t m, const size_t n)
    : state_space(n)
  {}

  /*
  template<>
  inline ss_vect<std::complex<double>>::ss_vect
  (const std::shared_ptr<mad::desc> d, const ord_t m, const size_t n)
    : state_space(n) {}
  */

  /*
    template<>
    inline ss_vect<tpsa>::ss_vect(const tpsa &t,  const size_t n)
    : state_space(n, tpsa(t, mad_tpsa_same)) {}
  */

  /*
  template<>
  inline ss_vect<tpsa>::ss_vect
  (const std::shared_ptr<mad::desc> d, const ord_t m, const size_t n)
    : state_space(n) {}
  */
   /*
   {
      // review this interface !
      // better just to clone
      // not all tpsa need to address the same number of variables or knobs

      this->state_space.reserve(n);
      for(size_t i=0; i<n; ++i) {
	this->state_space.push_back(std::move(gtpsa::tpsa(t, mad_tpsa_same)));
      }
    }
  */

  template<>
  inline void ss_vect<double>::cst(std::vector<double> &r) const {
    for (size_t i = 0; i < this->size(); ++i) { r[i] = this->state_space[i]; }
    // std::transform(begin(), this->state_space.end(), r.state_space.begin(),
    // 		   [](T& elem) { return elem.cst();});
  }

  template<>
  inline void ss_vect<double>::setConstantPartWithoutCheck
  (const std::vector<double> &r) {
    for (size_t i = 0; i < this->size(); ++i) { this->state_space[i] = r[i]; }
  }

  template<>
  inline void ss_vect<tpsa>::setConstantPartWithoutCheck
  (const std::vector<double> &r) {
    for (size_t i = 0; i < this->size(); ++i)
      { this->state_space[i].set(0, r[i]); }
  }

  /*
  template<>
  inline void ss_vect<std::complex<double>>::cst
  (std::vector<std::complex<double>>& r) const {
    for(size_t i = 0; i < this->size(); ++i) { r[i] = this->state_space[i]; }
    // std::transform
    //   (begin(), this->state_space.end(), r.state_space.begin(), [](T& elem)
    //   { return elem.cst();});
  }
  */

  template<>
  inline void ss_vect<tpsa>::set_zero(void) {
    std::for_each
      (this->state_space.begin(), this->state_space.end(), [](tpsa &v)
      { v.clear(); });
  }

  template<>
  inline void ss_vect<double>::set_zero(void) {
    for (auto &v: this->state_space) { v = 0.0; };
  }

  template<>
  inline void ss_vect<double>::_copyInPlace(const ss_vect<double> &o) {
    std::transform
      (o.state_space.begin(), o.state_space.end(), this->state_space.begin(),
       [](const double &elem) -> double { return elem; });
  }

  template<>
  inline ss_vect<double> ss_vect<double>::clone(void) const {
    this->checkSize(this->state_space);
    size_t size = this->state_space.size();
    std::vector<double> nvec;
    nvec.reserve(size);
    std::transform
      (this->state_space.begin(), this->state_space.end(), nvec.begin(),
       [](const double &o) {
      double r = o;
      return r;
    });
    ss_vect<double> nv =
      ss_vect(nvec[0], nvec[1], nvec[2], nvec[3], nvec[4], nvec[5]);
    return nv;
  }

  inline ss_vect<gtpsa::tpsa> operator+
  (const ss_vect<gtpsa::tpsa> &v1, const ss_vect<double> &v2) {
    auto r = v1.clone();
    r += v2;
    return r;
  }

  inline ss_vect<gtpsa::tpsa> operator-
  (const ss_vect<gtpsa::tpsa> &v1, const ss_vect<double> &v2) {
    auto r = v1.clone();
    r -= v2;
    return r;
  }

  /**
   * @todo: review if not implemented using setJacobian and
   * arma::mat(6, 6, arma::fill::eye)
   *
   */
  template<>
  inline void ss_vect<tpsa>::set_identity(void) {

    this->checkSize(this->state_space);
    const size_t n = this->state_space.size();

    for (size_t i = 0; i < n; ++i) {
      auto &t_tpsa = this->state_space[i];
      t_tpsa.clear();
      /* gtpsa starts counting variables at 1 ...*/
      t_tpsa.setVariable(0e0, i+1);
      // t_tpsa.setsm(std::vector<idx_t>{int(i+1), 1}, 0e0, 1e0);
    }
    // throw std::runtime_error("gtpsa set identity needs to be implemented!");
    // std::for_each(this->state_space.begin(), this->state_space.end(),
    // 		  [](gtpsa::tpsa& v) { v.clear(); };
  }

  /**
   * @brief start of derivative
   *
   * @param nv number of variables
   *
   * estimate the first argument for getv / setv
   * @todo recheck if not contained in Laurent's gtpsa,
   * @todo how to handle knob variables
   */
  size_t estimate_start_order(size_t order, size_t nv);

  /**
   * @brief collect all first deriviatives in a matrix
   *
   * @todo is that necessarily the jacobian
   *
   * It depends how the
   */
  template<>
  arma::mat ss_vect<tpsa>::jacobian(void) const;
  template<>
  void ss_vect<tpsa>::setJacobian(arma::mat &jac);

  template<>
  arma::cube gtpsa::ss_vect<tpsa>::hessian() const;
  template<>
  void gtpsa::ss_vect<tpsa>::setHessian(arma::cube& jac);


  /**
   * @brief first derivatives and cst term
   *
   * constant term at last session
   * @warning review if
   */
  template<>
  inline arma::mat ss_vect<tpsa>::toMatrix(void) {

    arma::mat mat(this->size(), this->size());

    //mat.fill(NAN);
    mat.fill(0e0);

    for (unsigned int j = 0; j < this->size(); j++) {
      auto &t_tpsa = this->state_space[j];
      for (unsigned int k = 0; k < this->size(); k++) {
	auto tmp = std::vector<idx_t>{int(k + 1), 1};
	if (t_tpsa.getDescription()->isvalidsm(tmp)) {
	  mat(j, k) = t_tpsa.getsm(tmp);
	} else {
	  std::cerr
	    << "toMatrix invalid k " << tmp[0] << ", o " << tmp[1] << "\n";
	}
      }
    }
    // mat(j, tps_n-1) = get_m_ij(map, j+1, 0);
    // mat(tps_n-1, tps_n-1) = 1e0;
    return mat;

  }

  /*
   * need to duplicate for ctpsa ...
   *
   * need to rearange ss_vect
   */
  template<>
  int ss_vect<tpsa>::getMaximumOrder(void) const;

  /**
   * @brief
   *
   * @return norm of tpsa function
   *
   * @todo: describe method
   */
  template<>
  double ss_vect<tpsa>::computeNorm(void) const;

  /* not inlined as quite some functionality behind the scenes */
  template<>
  void ss_vect<tpsa>::rcompose(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b);

  template<>
  void ss_vect<tpsa>::rminv(const ss_vect<tpsa>& a);

  template<>
  void ss_vect<tpsa>::rpminv
  (const ss_vect<tpsa>& a, std::vector<idx_t>& select);


  /**
   * @brief truncated power series to vector flow
   *
   * @todo: add derivative to flow
   * @param o: trunctated power series containing coefficients
   *
   *
   */
  template<>
  void ss_vect<tpsa>::rvec2fld(const tpsa& a);

  /**
   * @brief vector flow to truncated power series
   *
   * @param r: stores the computed truncated power serise in this object
   *
   * Intd in Forest's F77 LieLib.
   */

  template<>
  void ss_vect<tpsa>::fld2vec(tpsa * r) const;
  template<>
  void ss_vect<tpsa>::fgrad(tpsa * b, tpsa * r) const;
  template<>
  void ss_vect<tpsa>::rliebra(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b);
  template<>
  void ss_vect<tpsa>::rexppb(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b);
  template<>
  void ss_vect<tpsa>::rlogpb(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b);

  template<>
  void ss_vect<tpsa>::rgetOrder(const ss_vect<tpsa>& a, const int order);

  template<>
  void ss_vect<tpsa>::rderiv(const ss_vect<tpsa>& a, const int order);


  template<typename T>
  void get_mns
  (const ss_vect<T> &x, const int no1, const int no2, ss_vect<T> &y)
  { x.get_mns(no1, no2, y); }

  template<typename T>
  void M_to_h(const ss_vect<T> &M, T &h)
  { M.M_to_h(h); }

  template<typename T>
  void M_to_h_DF(const ss_vect<T> &M, T &h)
  { M.M_to_h_DF(h); }

  template<typename T>
  void GoFix(const ss_vect<T> &M, ss_vect<T> &A_0)
  { M.GoFix(A_0); }

  template<typename T>
  void Map_Norm
  (const ss_vect<T> &M, ss_vect<T> &A_0, ss_vect<T> &A_1, ss_vect<T> &R, T &g,
   T &K)
  { M.Map_Norm(A_0, A_1, R, g, K); }


  inline ss_vect<tpsa> compose(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b)
  {
    auto c = a.allocateLikeMe();
    c.rcompose(a, b);
    return c;
  }

  inline ss_vect<tpsa> liebra(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b)
  {
    auto r = a.allocateLikeMe();
    r.rliebra(a, b);
    return r;
  }

  inline ss_vect<tpsa> exppb(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b)
  {
    auto r = a.allocateLikeMe();
    r.rexppb(a, b);
    return r;
  }

  inline ss_vect<tpsa> logpb(const ss_vect<tpsa>& a, const ss_vect<tpsa>& b)
  {
    auto r = a.allocateLikeMe();
    r.rlogpb(a, b);
    return r;
  }

  inline ss_vect<tpsa> minv(const ss_vect<tpsa>& a) {
    auto r = a.allocateLikeMe();
    r.rminv(a);
    return r;
  }

  inline ss_vect<tpsa> pminv
  (const ss_vect<tpsa>& a, std::vector<idx_t>& select) {
    auto r = a.allocateLikeMe();
    r.rpminv(a, select);
    return r;
  }

} /* namespace gtpsa */

#endif /* _GTPSA_SS_VECT_H_ */

/*
 * Local Variables:
 * mode: c++
 * End:
 */
