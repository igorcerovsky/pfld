#ifndef FACET_H_
#define FACET_H_


#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <atomic>

namespace pfld {

#define EPS	1.0e-12
#define PI 3.1415926535897932384626433832795
#define PI2 6.283185307179586476925286766559
#define GRAVCONST 6.6738480e-11


template <typename T> T sign(T val) {	return (T(0) < val) - (val < T(0));}
//template <typename T> T sign(T val) { return copysign(1.0, val); }




template <typename T>
class Point3D
{
public:
	virtual ~Point3D() {};
	Point3D() : x(0), y(0), z(0) { };
	Point3D(T x_, T y_, T z_) : x(x_), y(y_), z(z_) { };
	Point3D(const Point3D& pt) : x(pt.x), y(pt.y), z(pt.z) {}
	const Point3D& operator = (const Point3D& pt)  { x = pt.x; y = pt.y; z = pt.z; return *this; }
	Point3D(const Point3D&& pt)
	{ x = std::move(pt.x); y = std::move(pt.y); z = std::move(pt.z); }
	const Point3D& operator = (const Point3D&& pt) 
	{ x = std::move(pt.x); y = std::move(pt.y); z = std::move(pt.z); return *this; }

	bool operator==(const Point3D& pt) const { return (x == pt.x && y == pt.y && z == pt.z); }
	bool operator!=(const Point3D& pt) const { return!( *this==pt); }

	Point3D operator-(const Point3D& pt) const { return Point3D(x - pt.x, y - pt.y, z - pt.z); }
	static void Sub(const Point3D& pt1, const Point3D& pt2, Point3D& res) 
	{
		res.x = pt1.x - pt2.x; res.y = pt1.y - pt2.y; res.z = pt1.z - pt2.z;
	}
	Point3D operator+(const Point3D& pt) const { return Point3D(x + pt.x, y + pt.y, z + pt.z); }
	static void Add(const Point3D& pt1, const Point3D& pt2, Point3D& res)
	{
		res.x = pt1.x + pt2.x; res.y = pt1.y + pt2.y; res.z = pt1.z + pt2.z;
	}

	void operator-=(const Point3D& pt) { x -= pt.x; y -= pt.y; z -= pt.z; }
	void operator+=(const Point3D& pt) { x += pt.x; y += pt.y; z += pt.z; }

	// vector arithmetics
	// scalar multiplication of two vectors in 3D
	T operator *(const Point3D& pt) const { return (x*pt.x + y*pt.y + z*pt.z); }
	// vector * number
	Point3D operator *(T a) const { return Point3D(a*x,a*y,a*z); };

	// vector cross product
	Point3D operator /(const Point3D& pt) const
	{
		return Point3D(y*pt.z - z*pt.y,
			z*pt.x - x*pt.z,
			x*pt.y - y*pt.x);
	}
	static void Cross(const Point3D& pt1, const Point3D& pt2, Point3D& res)
	{
		res.x = pt1.y*pt2.z - pt1.z*pt2.y;
		res.y = pt1.z*pt2.x - pt1.x*pt2.z;
		res.z = pt1.x*pt2.y - pt1.y*pt2.x;
	}
	// vector / number
	Point3D operator /(const T a) { return CPoint3D(x / a, y / a, z / a); }

	// vector length
	double Abs() const { return sqrt(x*x + y*y + z*z); };
	// unit vector
	void Unit() {
		double l = Abs();
		if (l > DBL_EPSILON) { x /= l; y /= l; z /= l; }
		else { x = 0; y = 0; z = 0; }
	}

	bool IsZero()
	{
		return (x < DBL_EPSILON && y < DBL_EPSILON && z < DBL_EPSILON);
	}

	bool IsEqualEps(const Point3D& pt, const double eps = DBL_EPSILON)
	{
		return (fabs(x - pt.x) < eps
			&& fabs(y - pt.y) < eps
			&& fabs(z - pt.z) < eps);
	}

	static bool AlmostEqual(const Point3D& pt1, const Point3D& pt2, const int digits=4)
	{
		return almost_equal(pt1.y, pt2.y, digits);
			//&& almost_equal(pt1.y, pt2.y, digits) 
			//&& almost_equal(pt1.z, pt2.z, digits);
	}

public:
	T x;
	T y;
	T z;
};



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