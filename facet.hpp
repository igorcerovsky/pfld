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


template <typename T> T sign(T val) {	return (T(0) < val) - (val < T(0));}
//template <typename T> T sign(T val) { return copysign(1.0, val); }


using valvec = std::vector < double >;
using point = Point3D < double >;
using ptvec = std::vector < Point3D<double> > ;
#ifdef FIELD_ATOMIC_DOUBLE
	using double_pfld = std::atomic < double > ;
#else
	using double_pfld = double;
#endif

class Facet
{
	enum class FacetType { normal, oposite };

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
	void Init(ptvec& pts, double densityCCW, double densityCW = 0.0);

	void Fld_G(const point &r, point& grv);
	void Fld_Gz(const point &r, double_pfld& g);
	void Fld_G(const point& r, const point& ro, const double& ro0, point& grv);
	void Fld_Gz(const point& r, const point& ro, const double& ro0, double& gz);

	ptvec& Data() { return _pts;}

	void FldGS(const point& r, const point M, point& mag, point& grv);
	void FldGS_M(const point& r, const point M, point& mag);
	void FldGS_G(const point& r, point& grv);
	void FldGS_Gz(const point& r, double& gz);
	double SolidAngle(ptvec& spts);

protected:
	void FldVlado(const point &r, double& f);
	void FldGS(const point& r, point& f);

};

using facet_vec = std::vector < Facet >;


}; // namespace pfld

#endif  // FACET_H_