#ifndef _GTPSA_PYTHON_OBJECTS_WITH_NAMED_INDEX_H_
#define _GTPSA_PYTHON_OBJECTS_WITH_NAMED_INDEX_H_ 1

#include <gtpsa/ctpsa.hpp>
#include <gtpsa/tpsa.hpp>
#include <gtpsa/ss_vect.h>
#include <gtpsa/python/name_index.h>
//#include "named_index.h"


namespace gtpsa::python {
  class ObjectWithNamedIndex {
    std::shared_ptr<gtpsa::python::IndexMapping> m_mapping;

  public:
    inline ObjectWithNamedIndex
    (std::shared_ptr<gtpsa::python::IndexMapping> mapping)
      : m_mapping(mapping)
    {}

    inline auto getMapping(void) const {
      return this->m_mapping;
    }
    inline void setMapping
    (std::shared_ptr<gtpsa::python::IndexMapping> mapping)
    {
      this->m_mapping = mapping;
    }
  };

  class TpsaWithNamedIndex : public gtpsa::tpsa, public ObjectWithNamedIndex {


  public:
    using base = gtpsa::tpsa;

    inline TpsaWithNamedIndex
    (std::shared_ptr<mad::desc> desc, const ord_t mo,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr): base(desc, mo)
      , ObjectWithNamedIndex(mapping)
    {}

    inline TpsaWithNamedIndex
    (const base& t,  std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr)
      : base(t)
      , ObjectWithNamedIndex(mapping)
    {}


    /* not accepting solely base object ... if mapping is lost, it is lost ...*/

    inline TpsaWithNamedIndex  clone(void) const
    { return TpsaWithNamedIndex( base::clone() ); }

    inline TpsaWithNamedIndex  operator  - ( void         ) const
    { return TpsaWithNamedIndex( base::operator-(), this->getMapping() ); }

    inline TpsaWithNamedIndex& operator += (const TpsaWithNamedIndex& o )
    { base::operator += (o) ; return *this; }
    inline TpsaWithNamedIndex& operator -= (const TpsaWithNamedIndex& o )
    { base::operator -= (o) ; return *this; }
    inline TpsaWithNamedIndex& operator *= (const TpsaWithNamedIndex& o )
    { base::operator *= (o) ; return *this; }
    inline TpsaWithNamedIndex& operator /= (const TpsaWithNamedIndex& o )
    { base::operator /= (o) ; return *this; }

    inline TpsaWithNamedIndex& operator += (const double o )
    { base::operator += (o) ; return *this; }
    inline TpsaWithNamedIndex& operator -= (const double o )
    { base::operator -= (o) ; return *this; }
    inline TpsaWithNamedIndex& operator *= (const double o )
    { base::operator *= (o) ; return *this; }
    inline TpsaWithNamedIndex& operator /= (const double o )
    { base::operator /= (o) ; return *this; }


    inline TpsaWithNamedIndex  operator +  ( const double o ) const
    { return TpsaWithNamedIndex( base::operator+ (o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator -  ( const double o ) const
    { return TpsaWithNamedIndex( base::operator- (o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator *  ( const double o ) const
    { return TpsaWithNamedIndex( base::operator* (o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator /  ( const double o ) const
    { return TpsaWithNamedIndex( base::operator/ (o), this->getMapping()  ) ; }

    inline TpsaWithNamedIndex  operator+ ( const TpsaWithNamedIndex& o ) const
    { return TpsaWithNamedIndex
	( static_cast<const gtpsa::tpsa&>(*this) +
	  static_cast<const gtpsa::tpsa&>(o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator- ( const TpsaWithNamedIndex& o ) const
    { return TpsaWithNamedIndex
	( static_cast<const gtpsa::tpsa&>(*this)
	  - static_cast<const gtpsa::tpsa&>(o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator* ( const TpsaWithNamedIndex& o ) const
    { return TpsaWithNamedIndex
	( static_cast<const gtpsa::tpsa&>(*this)
	  * static_cast<const gtpsa::tpsa&>(o), this->getMapping()  ) ; }
    inline TpsaWithNamedIndex  operator/ ( const TpsaWithNamedIndex& o ) const
    { return TpsaWithNamedIndex
	( static_cast<const gtpsa::tpsa&>(*this)
	  / static_cast<const gtpsa::tpsa&>(o), this->getMapping()  ) ; }

#if 0
    friend inline TpsaWithNamedIndex  operator +
    ( const TpsaWithNamedIndex&  a,  const TpsaWithNamedIndex&  b );
    friend inline TpsaWithNamedIndex  operator -
    ( const TpsaWithNamedIndex&  a,  const TpsaWithNamedIndex&  b );
    friend inline TpsaWithNamedIndex  operator *
    ( const TpsaWithNamedIndex&  a,  const TpsaWithNamedIndex&  b );
    friend inline TpsaWithNamedIndex  operator /
    ( const TpsaWithNamedIndex&  a,  const TpsaWithNamedIndex&  b );
#endif
  };

  inline TpsaWithNamedIndex  operator +
  ( const double a,  const TpsaWithNamedIndex&  b)
  { return TpsaWithNamedIndex( a + gtpsa::tpsa(b), b.getMapping()  ) ; }
  inline TpsaWithNamedIndex  operator -
  ( const double a,  const TpsaWithNamedIndex&  b)
  { return TpsaWithNamedIndex( a - gtpsa::tpsa(b), b.getMapping()  ) ; }
  inline TpsaWithNamedIndex  operator *
  ( const double a,  const TpsaWithNamedIndex&  b)
  { return TpsaWithNamedIndex( a * gtpsa::tpsa(b), b.getMapping()  ) ; }
  inline TpsaWithNamedIndex  operator /
  ( const double a,  const TpsaWithNamedIndex&  b)
  { return TpsaWithNamedIndex( a / gtpsa::tpsa(b), b.getMapping()  ) ; }

  inline TpsaWithNamedIndex pow
  (const TpsaWithNamedIndex& a,  const TpsaWithNamedIndex& b)
  { return TpsaWithNamedIndex
      ( gtpsa::pow
	(static_cast<const TpsaWithNamedIndex::base&>(a),
	 static_cast<const TpsaWithNamedIndex::base&>(b) ) ); }
  inline TpsaWithNamedIndex pow
  (const TpsaWithNamedIndex& a,  const int   i)
  { return TpsaWithNamedIndex
      ( gtpsa::pow(static_cast<const TpsaWithNamedIndex::base&>(a), i) ); }
  inline TpsaWithNamedIndex pow
  (const TpsaWithNamedIndex& a,  const num_t v)
  { return TpsaWithNamedIndex
      ( gtpsa::pow(static_cast<const TpsaWithNamedIndex::base&>(a), v  )); }


#define GTPSA_FUNC_ARG1_WITH_RET_ARG(fname)				\
  inline void fname ## _ (const TpsaWithNamedIndex& t, TpsaWithNamedIndex* r) \
  { r->apply_with_return_object(t, mad::fname);  }
#define GTPSA_FUNC_ARG1_NO_RET_ARG(fname)				\
  inline TpsaWithNamedIndex fname (const TpsaWithNamedIndex& t)		\
  { return apply<TpsaWithNamedIndex>(t, fname ## _);  }
#define GTPSA_FUNC_ARG1(fname) GTPSA_FUNC_ARG1_WITH_RET_ARG(fname)	\
  GTPSA_FUNC_ARG1_NO_RET_ARG(fname)

#include <gtpsa/funcs.h>
#undef GTPSA_FUNC_ARG1_NO_RET_ARG
#undef GTPSA_FUNC_ARG1_WITH_RET_ARG
#undef GTPSA_FUNC_ARG1


  class CTpsaWithNamedIndex : public gtpsa::ctpsa {
    std::shared_ptr<gtpsa::python::IndexMapping> m_mapping;

  public:
    using base = gtpsa::ctpsa;

    CTpsaWithNamedIndex
    (std::shared_ptr<mad::desc> desc, const ord_t mo,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(desc, mo), m_mapping(mapping)
    {}

    CTpsaWithNamedIndex
    (const tpsa& t, const ord_t mo,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(t, mo), m_mapping(mapping)
    {}

    CTpsaWithNamedIndex
    (const base& t,  std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(t), m_mapping(mapping)
    {}

    inline auto getMapping(void) const {
      return this->m_mapping;
    }
    inline void setMapping
    (std::shared_ptr<gtpsa::python::IndexMapping> mapping) {
      this->m_mapping = mapping;
    }

    /* not accepting solely base object ... if mapping is lost, it is lost ...*/

    inline CTpsaWithNamedIndex  clone(void) const
    { return CTpsaWithNamedIndex( base::clone(), this->getMapping()  ); }
    inline CTpsaWithNamedIndex  newFromThis(void) const
    { return CTpsaWithNamedIndex( base::newFromThis(), this->getMapping()  ); }

    inline CTpsaWithNamedIndex  operator  - ( void         ) const
    { return CTpsaWithNamedIndex( base::operator-() ); }

    inline CTpsaWithNamedIndex& operator += (const CTpsaWithNamedIndex& o )
    { base::operator += (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator -= (const CTpsaWithNamedIndex& o )
    { base::operator -= (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator *= (const CTpsaWithNamedIndex& o )
    { base::operator *= (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator /= (const CTpsaWithNamedIndex& o )
    { base::operator /= (o) ; return *this; }

    inline CTpsaWithNamedIndex& operator += (const std::complex<double> o )
    { base::operator += (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator -= (const std::complex<double> o )
    { base::operator -= (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator *= (const std::complex<double> o )
    { base::operator *= (o) ; return *this; }
    inline CTpsaWithNamedIndex& operator /= (const std::complex<double> o )
    { base::operator /= (o) ; return *this; }


    inline CTpsaWithNamedIndex  operator +  ( const std::complex<double> o )
      const
    { return CTpsaWithNamedIndex( base::operator+ (o), this->getMapping()  ) ;
    }
    inline CTpsaWithNamedIndex  operator -  ( const std::complex<double> o )
      const
    { return CTpsaWithNamedIndex( base::operator- (o), this->getMapping()  ) ;
    }
    inline CTpsaWithNamedIndex  operator *  ( const std::complex<double> o )
      const
    { return CTpsaWithNamedIndex( base::operator* (o), this->getMapping()  ) ; }
    inline CTpsaWithNamedIndex  operator /  ( const std::complex<double> o )
      const
    { return CTpsaWithNamedIndex( base::operator/ (o), this->getMapping()  ) ; }

  };

  inline CTpsaWithNamedIndex  operator+
  ( const CTpsaWithNamedIndex& a, const CTpsaWithNamedIndex& b )
  { return CTpsaWithNamedIndex
      ( static_cast<const gtpsa::ctpsa&>(a)
	+ static_cast<const gtpsa::ctpsa&>(b), a.getMapping()  ) ; }
  inline CTpsaWithNamedIndex  operator-
  ( const CTpsaWithNamedIndex& a, const CTpsaWithNamedIndex& b )
  { return CTpsaWithNamedIndex
      ( static_cast<const gtpsa::ctpsa&>(a)
	- static_cast<const gtpsa::ctpsa&>(b), a.getMapping()  ) ; }
  inline CTpsaWithNamedIndex  operator*
  ( const CTpsaWithNamedIndex& a, const CTpsaWithNamedIndex& b )
  { return CTpsaWithNamedIndex
      ( static_cast<const gtpsa::ctpsa&>(a)
	* static_cast<const gtpsa::ctpsa&>(b), a.getMapping()  ) ; }
  inline CTpsaWithNamedIndex  operator/
  ( const CTpsaWithNamedIndex& a, const CTpsaWithNamedIndex& b )
  { return CTpsaWithNamedIndex
      ( static_cast<const gtpsa::ctpsa&>(a)
	/ static_cast<const gtpsa::ctpsa&>(b), a.getMapping()  ) ; }

  inline CTpsaWithNamedIndex pow
  (const CTpsaWithNamedIndex& a,  const CTpsaWithNamedIndex& b)
  { return CTpsaWithNamedIndex
      ( gtpsa::pow(static_cast<const CTpsaWithNamedIndex::base&>(a),
		   static_cast<const CTpsaWithNamedIndex::base&>(b) ) ); }
  inline CTpsaWithNamedIndex pow
  (const CTpsaWithNamedIndex& a,  const int   i)
  { return CTpsaWithNamedIndex
      ( gtpsa::pow(static_cast<const CTpsaWithNamedIndex::base&>(a), i) ); }
  inline CTpsaWithNamedIndex pow
  (const CTpsaWithNamedIndex& a,  const num_t v)
  { return CTpsaWithNamedIndex
      ( gtpsa::pow(static_cast<const CTpsaWithNamedIndex::base&>(a), v  )); }


#define GTPSA_FUNC_ARG1_WITH_RET_ARG(fname)				\
  inline void fname ## _						\
  (const CTpsaWithNamedIndex& t, CTpsaWithNamedIndex* r)		\
  { r->apply_with_return_object(t, mad::fname);  }
#define GTPSA_FUNC_ARG1_NO_RET_ARG(fname)				\
  inline CTpsaWithNamedIndex fname (const CTpsaWithNamedIndex& t)	\
  { return apply<CTpsaWithNamedIndex>(t, fname ## _);  }
#define GTPSA_FUNC_ARG1(fname) GTPSA_FUNC_ARG1_WITH_RET_ARG(fname)	\
  GTPSA_FUNC_ARG1_NO_RET_ARG(fname)

#include <gtpsa/funcs.h>
#undef GTPSA_FUNC_ARG1_NO_RET_ARG
#undef GTPSA_FUNC_ARG1_WITH_RET_ARG
#undef GTPSA_FUNC_ARG1


  /**
   * State space with indices that can be used as attributed
   *
   * Note:
   *    Enable shared is required as the python wrappers need to provide .loc
   *    and .iloc
   *    These return new objects that need to take a shared pointer to this
   *    object.
   *    See :meth:`getPtr`
   */

  template<typename T>
  class StateSpaceWithNamedIndex :
    public gtpsa::ss_vect<T>,
    public std::enable_shared_from_this<StateSpaceWithNamedIndex<T>> {
    std::shared_ptr<gtpsa::python::IndexMapping> m_mapping;

    using base = typename gtpsa::ss_vect<T>;

  public:
    ~StateSpaceWithNamedIndex() {
      // std::cerr << "deleting StateSpaceWithNamedIndex: " << this
      // 		<< std::endl;
      // std::cerr.flush();
    }
    StateSpaceWithNamedIndex
    (const std::shared_ptr<gtpsa::mad::desc> desc, const ord_t mo,
     const size_t n = ss_vect_n_dim,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(desc, mo, n), m_mapping(mapping)
    { /* std::cerr << "allocating StateSpaceWithNamedIndex: " << this
		<< std::endl; */ }

    StateSpaceWithNamedIndex
    (const T& t, const size_t n = ss_vect_n_dim,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(t, n), m_mapping(mapping)
    { /* std::cerr << "allocating StateSpaceWithNamedIndex: " << this
		<< std::endl; */ }


    StateSpaceWithNamedIndex
    (const base& vec, std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) : base(vec), m_mapping(mapping)
    {}

    /*
    StateSpaceWithNamedIndex
    (const std::vector<T> &vec, const size_t n = ss_vect_n_dim,
     std::shared_ptr<gtpsa::python::IndexMapping> mapping =
     gtpsa::python::default_index_mapping_ptr) :
      base(vec, n), m_mapping(mapping)
    {}
    */


    inline auto getPtr(void) {
      return std::enable_shared_from_this<StateSpaceWithNamedIndex<T>>::
	shared_from_this();
    }
    inline auto getMapping(void) const {
      return this->m_mapping;
    }
    inline void setMapping(std::shared_ptr<gtpsa::python::IndexMapping> mapping)
    {
      this->m_mapping = mapping;
    }

    inline StateSpaceWithNamedIndex clone(void) const
    { return StateSpaceWithNamedIndex(base::clone(), this->getMapping()); }
    inline StateSpaceWithNamedIndex<double> cst(void) const
    { return StateSpaceWithNamedIndex<double>
	(base::cst(), this->getMapping()); }

    inline StateSpaceWithNamedIndex& operator +=
    (const StateSpaceWithNamedIndex& o) { base::operator += (o); return *this; }
    inline StateSpaceWithNamedIndex& operator -=
    (const StateSpaceWithNamedIndex& o) { base::operator -= (o); return *this; }

    template<typename U>
    inline StateSpaceWithNamedIndex<T>& operator +=
    (const StateSpaceWithNamedIndex<U>& o)
    { base::operator += (o); return *this; }
    template<typename U>
    inline StateSpaceWithNamedIndex<T>& operator -=
    (const StateSpaceWithNamedIndex<U>& o)
    { base::operator -= (o); return *this; }

    inline StateSpaceWithNamedIndex& operator +=
    (const double o) { base::operator +=(o); return *this; }
    inline StateSpaceWithNamedIndex& operator -=
    (const double o) { base::operator -=(o); return *this; }

    inline StateSpaceWithNamedIndex  operator - (void) const
    { return StateSpaceWithNamedIndex
	( base::operator- (), this->getMapping()); }
    inline StateSpaceWithNamedIndex  operator +
    (const StateSpaceWithNamedIndex& o) const
    { return StateSpaceWithNamedIndex
	( base::operator+ (o), this->getMapping()); }
    inline StateSpaceWithNamedIndex  operator -
    (const StateSpaceWithNamedIndex& o) const
    { return StateSpaceWithNamedIndex
	( base::operator- (o), this->getMapping()); }

    // inline StateSpaceWithNamedIndex operator *
    // (const StateSpaceWithNamedIndex& o) const
    // { return StateSpaceWithNamedIndex
    // 	(base(*this) * base(o), this->getMapping());}
    // inline StateSpaceWithNamedIndex operator /
    // (const StateSpaceWithNamedIndex& o) const
    // { return StateSpaceWithNamedIndex
    // 	(base(*this) / base(o), this->getMapping());}

    inline StateSpaceWithNamedIndex  operator +
    (const double o) const
    { return StateSpaceWithNamedIndex
	(base::operator+ (o), this->getMapping()); }
    inline StateSpaceWithNamedIndex  operator - (const double o) const
    { return StateSpaceWithNamedIndex
	(base::operator- (o), this->getMapping()); }
  };

  /*
    still required ?
  inline StateSpaceWithNamedIndex<gtpsa::tpsa>& operator +=
  (StateSpaceWithNamedIndex<gtpsa::tpsa>& a,
   const StateSpaceWithNamedIndex<double>& b)
  {(ss_vect<gtpsa::tpsa>&)(a) += (b);	return a;  }
  inline StateSpaceWithNamedIndex<gtpsa::tpsa>& operator -=
  (StateSpaceWithNamedIndex<gtpsa::tpsa>& a,
   const StateSpaceWithNamedIndex<double>& b)
  {(ss_vect<gtpsa::tpsa>&)(a) -= (b);	return a;  }
  */

  inline StateSpaceWithNamedIndex<gtpsa::tpsa> operator+
  (const StateSpaceWithNamedIndex<gtpsa::tpsa> &v1,
   const StateSpaceWithNamedIndex<double> &v2) {
    auto r = v1.clone();
    r += v2;
    return r;
  }

  inline StateSpaceWithNamedIndex<gtpsa::tpsa> operator-
  (const StateSpaceWithNamedIndex<gtpsa::tpsa> &v1,
   const StateSpaceWithNamedIndex<double> &v2) {
    auto r = v1.clone();
    r -= v2;
    return r;
  }

} // namespace gtpsa::python

#endif /* _GTPSA_PYTHON_OBJECTS_WITH_NAMED_INDEX_H_ */
