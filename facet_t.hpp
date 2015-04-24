#ifndef FACET_H_
#define FACET_H_


#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <atomic>

#include "point.hpp"


namespace pfld {

#define EPS	1.0e-12
#define PI 3.1415926535897932384626433832795
#define PI2 6.283185307179586476925286766559
#define GRAVCONST 6.6738480e-11


template <typename T> T sign(T val) {	return T(0. < val) - (val < 0.);}
//template <typename T> T sign(T val) { return copysign(1.0, val); }



template<typename T>
class Facet
{
public:
	using valvec = std::vector < T >;
	using point = Point3D < T >;
	using ptvec = std::vector < Point3D<T> >;

#ifdef FIELD_ATOMIC_DOUBLE
	using double_pfld = std::atomic < T >;
#else
	using double_pfld = T;
#endif

public:
	virtual ~Facet() {};
	Facet();
	Facet(const ptvec& pts);
	Facet(const Facet& fct);

	Facet& operator=(const Facet& fct);
	bool operator==(const Facet& fct) const;

	void operator()(const point& r, point& grv);
	void operator()(const point& r, double_pfld& g);

protected:
	bool   _initialized; // is initialized?
	size_t _sz;
	int    _id;    // facet ID
	ptvec  _pts;   // points
	ptvec  _L;     // array of vector l vectors; l = pts[i+1] - pts[i]	
	ptvec  _mi;    // array of unit vector mi vectors; mi = Unit(pts[i+1] - pts[i])	
	ptvec  _ni;    // arrys of unit ni vectors; ni = _mi[i] / _n;
	point  _n;     // facet unit normal vector
	valvec _len;   // arry of side lengths


public:
	void Init();
	void Init(ptvec& pts);

	void Fld_G(const point &r, point& grv);
	void Fld_Gz(const point &r, double_pfld& g);
	void Fld_G(const point& r, const point& ro, const T& ro0, point& grv);
	void Fld_Gz(const point& r, const point& ro, const T& ro0, T& gz);

	void FldGS(const point& r, const point M, point& mag, point& grv);
	void FldGS_M(const point& r, const point M, point& mag);
	void FldGS_G(const point& r, point& grv);
	void FldGS_Gz(const point& r, T& gz);

	static T SolidAngle(ptvec& pts, const T inOut, const size_t sz);

protected:
	void FldVlado(const point &r, T& f);
	void FldGS(const point& r, point& f);

};



}; // namespace pfld

#endif  // FACET_H_