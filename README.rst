C++ Interface to Truncated Power Series Algebra Module Interface
=================================================================

This code base is providing a shallow c++ wrapper to the
truncatad power series algebra module as provided in mad-ng

	https://github.com/MethodicalAcceleratorDesign/MAD.

For details see

	https://github.com/MethodicalAcceleratorDesign/MAD/blob/dev/src/libgtpsa/README.GTPSA.

	http://mad.web.cern.ch/mad/releases/madng/html/index.html

	http://mad.web.cern.ch/mad/releases/madng/html/mad_mod_diffalg.html


**NB**: this code base modifies the C function defintions of the original code.
For serious work please checkout the original code, in particular if you are using the "C" language.

C++ and Python Pybind11 Interfaces
==================================

The Python Pybind11 <- C++ <- C *gtpsa* interface was prototyped and implemented by Pierre Schnizer.

References:
	P\. Schnizer, W. Khail, J. Bengtsson *Small Talk on AT* IPAC 2022

	https://accelconf.web.cern.ch/ipac2022/papers/tupost029.pdf

	L\. Deniau, C. Tomoiagă *Generalised Truncated Power Series Algebra for Fast Particle Accelerator
	Transport Maps* IPAC 2015

	https://accelconf.web.cern.ch/ipac2015/papers/mopje039.pdf

Turned out that the CEERN gtpsa map concatenator can not handle parameter dependence; so it had to be
reimplemented.

The C++ <- C gtpsa bridge interface is in:

	../src/gtpsa/python/src/gtpsa.cc

	/*
	 * Implementation split up in three parts:
	 *
	 * 1. bridge: defined as template
	 * 2. operator functions using the bridge
	 * 3. class using c++ operators defined in a template
	 * 4. class providing full functionality derived from the template
	 *
	 * This splits the functionality in different parts. Hopefully that
	 * makes the code a little more maintainable
	 */


However, some of the key *gtpsa* map analysis functions are implemented in the *Lua* scripting language;
see below.

Hence, they have been re-implemented in C++.

Data Types
----------

	| num_t double
	| ord_t unsigned char
	| idx_t int32_t
	| ssz_t int32_t



Interfacing gtpsa Functions
---------------------------

E.g. deriv.

Files:

	../src/gtpsa/c++/gtpsa/tpsa.hpp
	../src/gtpsa/c++/gtpsa/bridge/brigde.hpp
	../src/gtpsa/c++/gtpsa/intern/with_operators.h
	../src/gtpsa/python/src/gtpsa_delegator.h
	../src/gtpsa/python/src/gtpsa.cc

C++ -> Python Pybind11 Part
---------------------------
The *gtpsa* Python Pybind11 <- C++ part is in:

Python interface:

	../python/src/thor_scsi.cc

		| M_to_h_DF
		| h_DF_to_M
		| CtoR
		| RtoC
		| GoFix
		| Map_Norm

and the implementation:

	../src/gtpsa/c++/gtpsa/lielib.cc

	../src/gtpsa/python/src/desc.cc

		| number_of_variables(ord_t \*mo_=0, int \*np_=0, ord_t \*po_=0) -> int
		| maximum_orders(int nn=0, ord_t \*no=nullptr) -> int
		| maximum_length(ord_t mo) -> int
		| mono(idx_t i, std::vector<ord_t> \*m) -> int
		| indexsm -> int
		| # E.g.:
		|   print(desc.number_of_variables(0, 0, 0))
		|
		|   exps = np.zeros(nv, dtype=int)
		|   ord = desc.mono(0, exps)
		|   print(ord, exps)
		|
		|   print(desc.index([1, 0, 0, 0, 0, 0, 0]))

	../src/gtpsa/python/src/ss_vect.cc

		| # Support a .loc["x"] access to the elements.
		|     template<class WrappedClass, class P_MGR, typename T>
		|
		| # print (__str__) calls:
		| pstr
		|
		| iloc[]
		| # E.g.:
		|     map.iloc[k]
		| getOrder
		| set_zero(void)
		| truncate
		| # E.g.:
		|     desc.truncate(3)

	TPSA map operations:

		| deriv
		| (integ)
		| mnrm
		| fld2vec
		| fgrad
		| liebra
		| exppb
		| logpb
		| compose
		| inv
		| pinv

	../src/gtpsa/python/src/gtpsa.cc
and
	../src/gtpsa/python/src/gtpsa_delegator.h

		| # For functions returning a tpsa.
		|
		| print
		| (Sets *eps* 1e-30 vs. 0 for the *gtpsa* print function to supress printing of zeroes)
		| length
		| get_description
		| # E.g.:
		|     print(a.get_description())
		| get
		| set
		| getv
		| setv

		...

The *gtpsa* C++ <- C functions are in:

	../src/gtpsa/c++/gtpsa/python/objects_with_named_index.h

		| Basis arithmetic operators: [+, -, *, /,...].

	../src/gtpsa/c++/gtpsa/bridge/bridge.hpp

		| clear(void)

		| # Parameters: (constant part, monomial index, value).
		| setVariable(const base_type v, const idx_t iv = 0, const base_type scale = 0).

		| # Return order & exponents for monomial with index i.
		| mono(idx_t i, std::vector<ord_t> \*m) -> int
		
		| # Return index for monomial m.
		| #   use string for the exponents:
		| index(std::string s)
		| #   use array for the exponents:
		| index(const std::vector<ord_t> &m)
		| #   sparse monomials [(i, o)]:
		| indexsm(const std::vector<int> m)
		
		| # Return a pair (.first, .second) for ???
		| cycle(const idx_t i, std::vector<ord_t> \*m)

		| cst()

		| # Get constant term.
		| get(void)                           get()

		| # Get mon
		| get(const idx_t i)                  get(46)
		| get(const std::string s)            get()
		| get(const std::vector<ord_t> &m)    get(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0})
		| getsm(const std::vector<int> &m)

		| set(void)
		| ...

		| # The 1st parameter is the offset - set to 1, to skip constant part: 0..
		| getv(idx_t i, std::vector<base_type> *v)
		| setv(idx_t i, const std::vector<base_type> &v)

		| rgetorder
		| print
		| cst
		| pow
		| add
		| dif
		| sub
		| mul
		| div
		| acc
		| scl
		| inv
		| invsqrt
		| rderiv
		| rinteg
		| ...

	../src/gtpsa/c++/gtpsa/mad/wrapper.tpp

		| print()
		| print("", 1e-30, 0, stdout) (For TPSA vector; use cout << for map)
		|
		| rgetOrder
		|
		| setvar(const GTPSA_BASE_T v, const idx_t iv = 0, const GTPSA_BASE_T scl = 0)
		| mono(const idx_t i, std::vector<ord_t> \*m)
		| idxs(const std::string s)
		| idxm(const std::vector<ord_t> &m)
		| idxsm(const std::vector<int> m)
		| cycle(const idx_t i, std::vector<ord_t> \*m, GTPSA_BASE_T \*v)
		|
		| get0(void)                           get()
		| geti(const idx_t i)                  get(46)
		| gets(const std::string s)            get()
		| getm(const std::vector<ord_t> &m)    get(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0})
		| getsm(const std::vector<int> &m)
		|
		| # The 1st parameter is offset - 1 to skip constant part: 0..
		| getv(const idx_t i, std::vector<GTPSA_BASE_T> \*v)
		| setv(const idx_t i, const std::vector<GTPSA_BASE_T> &v)
		|
		| a*x[0]+b
		| set0(const num_t a, const num_t b)
		|
		| a*x[i]+b
		| seti(const idx_t i, const num_t a, const num_t b)
		|
		| a*x[m]+b
		| sets(const std::string &s, const num_t a, const num_t b)
		|
		| a*x[m]+b
		| setm(const std::vector<ord_t> &m, const num_t a, const num_t b)
		|
		| rderiv
		| rinteg

	../src/gtpsa/c++/gtpsa/mad/tpsa_wrapper.hpp
	Wrapper for C++ <- C.

		| norm
		| equ

	../src/gtpsa/c++/gtpsa/bridge/container.hpp

		| size
		| getMaximumOrder
		| computeNorm
		| rvec2fld
		| ...

	../src/gtpsa/c++/gtpsa/mad/container_wrapper.tpp

		| size
		| getMaximumOrder
		| computeNorm
		| rvec2fld
		| fld2vec
		| fgrad
		| rliebra
		| rexppb
		| rlogpb
		| rcompose (which call compose in the gtpsa library)
		| rminv
		| rpminv

	../src/gtpsa/c++/gtpsa/intern/with_operators.hpp

		| # The Python interface for maps calls:
		| pstr
		| # which calls:
		| show()
		| # For TPSA vector: only prints leading order - *level* parameter not implemented.
		| show(stdout, level)
		| print("", eps, 0)
		| operator<<


The *gtpsa* print functions are in:

	../src/gtpsa/mad-ng/src/mad_tpsa.c
	
		| mad_tpsa_setvar(tpsa_t \*t, num_t v, idx_t iv, num_t scl)
		|
		| mad_tpsa_mono(const tpsa_t \*t, idx_t i,  ssz_t n, ord_t m[])
		| mad_tpsa_idxs(const tpsa_t \*t, ssz_t n, str_t s)
		| mad_tpsa_idxm(const tpsa_t \*t, ssz_t n, const ord_t m[])
		| mad_tpsa_idxsm(const tpsa_t \*t, ssz_t n, const int m[])
		| mad_tpsa_cycle(const tpsa_t \*t, idx_t i, ssz_t n, ord_t m[], num_t \*v)
		|
		| mad_tpsa_get0(const tpsa_t \*t)
		| mad_tpsa_geti(const tpsa_t \*t, idx_t i)
		| mad_tpsa_gets(const tpsa_t \*t, ssz_t n, str_t s)
		| mad_tpsa_getm(const tpsa_t \*t, ssz_t n, const ord_t m[])
		| mad_tpsa_getsm(const tpsa_t \*t, ssz_t n, const int m[])
		|
		| # The 2nd parameter is offset - 1 to skip constant part: 0..
		| mad_tpsa_getv(const tpsa_t \*t, idx_t i, ssz_t n, num_t v[])

	../src/gtpsa/mad-ng/src]/mad_tpsa_io.c

	../src/gtpsa/mad-ng/src]/mad_tpsa_comp.c

		| print
		| print_damap

*Gtpsa* C++ <- C Interface
--------------------------
The general *gtpsa* C++ <- C interface is in:

	../src/gtpsa/c++/gtpsa/desc.hpp

	../src/gtpsa/c++/gtpsa/desc.cc

		| show
		| # Prints out info, e.g.:
		| #   id=2, nn=7, nv=7, np=0, mo=5, po=0, to=5, uno=0, no=[5555555]
		| info(FILE * fp = nullptr)
		| 
		| getDescription()->
		|    # Get all the info:
		|      getInfo
		|    #  e.g.:
		|      .getDescription()->getInfo()
		|    getNv
		|    maxOrd
		|    maxLen
		|
		| getNumberOfVariables
		| getVariablesMaximumOrder
		| getNumberOfParameters
		| getParametersMaximumOrder
		| getTotalNumber
		| getOrderPerParameter
		| getNv(ord_t \*mo_=0, int \*np_=0, ord_t \*po_=0)
		| maxOrd(int nn=0, ord_t \*no=nullptr)
		| maxLen(ord_t mo)
		| # Sets *to* for all gtpsa elements???
		| trunc(const ord_t to)
		| # E.g.:
		|     .getDescription()->trunc(k)

	../src/gtpsa/c++/gtpsa/ss_vect.cc

		| # For functions returning an ss_vect<>.
		|
		| # For general indexing:
		|     idx()
		|
		| ss_vect_n_dim
		| ss_vect
		| state_space
		| # For TPSA map: only prints leading order - *level* parameter not implemented.
		| show(std::ostream &strm, int level = 1, bool with_endl = true)
		|
		| jacobian
		| hessian
		| set_zero
		| set_identity
		| setConstant
		| setJacobian
		| setHessian
		| rcompose

	../src/gtpsa/c++/gtpsa/funcs.h

		| sqrt
		| exp
		| log
		| ...

	../src/gtpsa/c++/gtpsa/lielib.cc

		| M_to_h_DF
		| h_DF_to_M
		| CtoR
		| RtoC
		| GoFix
		| Map_Norm


TPSA descriptor operations:

	../src/gtpsa/mad-ng/src/mad_desc.h

	../src/gtpsa/mad-ng/src/mad_desc.c

		| int mad_desc_getnv(const D \*d, ord_t \*mo_, int \*np_, ord_t \*po_)
		| ord_t mad_desc_maxord(const D \*d, int n, ord_t no_[n])
		| # Sets *to* for all gtpsa elements???
		| ord_t mad_desc_gtrunc(const desc_t \*d, ord_t to)
		| void mad_desc_info(const D \*d, FILE \*fp_)

TPSA vector operations:

	../src/gtpsa/mad-ng/src/mad_tpsa.h

	../src/gtpsa/mad-ng/src/mad_tpsa_ops.c

		| add
		| sub
		| ...
		| integ
		| deriv
		| poisbra
		| ...
		| print
		| ...
		| cutord

TPSA map operations:

	../src/gtpsa/mad-ng/src/mad_tpsa_comp.c

		| Local
		|
		| print_damap
		|
		| Public
		|
		| compose
		| translate
		| eval


	../src/gtpsa/mad-ng/src]/mad_tpsa_comp_s.tc

		| compose

	../src/gtpsa/mad-ng/src]/mad_tpsa_minv.c

		| minv
		|
		| pinv

	../src/gtpsa/mad-ng/src/mad_tpsa_mops.c

		| Local
		|
		| print_damap
		|
		| Public
		|
		| exppb
		| logpb
		| liebra
		| fgrad
		|
		| Compute (Eq. (34)):
			| G(x;0) = -J grad.f(x;0)
		| vec2fld
		|
		| Compute(Eqs. (34)-(37)):
			| f(x;0) = \int_0^x J G(x';0) dx' = x^t J phi G(x;0)
		|
		| fld2vec
		| mnrm (norm)

Also, a few are in:

(coded in *Lua*)

	../src/gtpsa/mad-ng/src/madl_damap.mad

		| map_ctor
		| factor_map
		|
		| Factored Lie of exponential and poisson bracket:
		|
			| r = exp(:y1:) exp(:y2:)... x
		|
		| lieexppb
		| flofacg
		| ...

	../src/gtpsa/madl_gphys.mad

		| make_symp (Make map symplectic, thesis by Liam Healy)
		|
			| L\. Healy *Lie-Algebraic Methods for Treating Lattice Parameter Errors in Particle
			| Accelerators* Thesis, Univ. of Maryland, 1986.
		|
		| gphys.normal_ng (Map normal form)
		| normal_c        (Phasor basis)

*Lua* Scripting Language
----------------------
The *Lua* scripting language (Portuguese: *lua* -> *moon*) was created by the Computer Graphics
Technology Group (Tecgraf) at the PUC Univ., Rio de Janeiro, Brazil in 1993:

	https://www.lua.org/about.html

LuaJiT is a just-in-time compiler:

	https://luajit.org/luajit.html
