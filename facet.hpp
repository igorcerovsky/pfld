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
#define GRAVCONST 6.6738480e-11


//template<typename T> T sign(T d) { return ((d == 0) ? 0.0 : ((d<0) ? -1 : 1)); }


template <typename T> T sign(T val) {	return (T(0) < val) - (val < T(0));}

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
	Point3D operator+(const Point3D& pt) const { return Point3D(x + pt.x, y + pt.y, z + pt.z); }

	Point3D operator-=(const Point3D& pt) { x -= pt.x; y -= pt.y; z -= pt.z; }
	void operator+=(const Point3D& pt)    { x += pt.x; y += pt.y; z += pt.z; }

	// vector arithmetics
	// scalar multiplication of two vectors in 3D
	T operator *(const Point3D& pt) const { return (x*pt.x + y*pt.y + z*pt.z); }
	// vector * number
	Point3D operator *(T a) const { return Point3D(a*x,a*y,a*z); };

	// vector multiplication of two vectors in 3D
	Point3D operator /(const Point3D& pt) const 
	{
		return Point3D(y*pt.z - z*pt.y,
			z*pt.x - x*pt.z,
			x*pt.y - y*pt.x);
	}
	// vector / number
	Point3D operator /(const T a) { return CPoint3D(x / a, y / a, z / a); }

	// another operator maybe useful sometimes
	Point3D operator+(const T a) const { return CPoint3D(x + a, y + a, z + a); }

	// vector length
	double Abs() const { return sqrt(x*x + y*y + z*z); };
	// unit vector
	void Unit() {
		auto eps = std::numeric_limits<T>::epsilon();
		double l = Abs();
		if (l > eps ) { x /= l; y /= l; z /= l; }
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

public:
	T x;
	T y;
	T z;
};


enum class FacetType { normal, oposite };

using valvec = std::vector < double >;
using point = Point3D < double >;
using ptvec = std::vector < Point3D<double> > ;
using atmic_double = std::atomic < double > ;

class Facet
{

public:
	virtual ~Facet() {};
	Facet();
	Facet(const ptvec& pts);
	Facet(const Facet& fct);

	Facet& operator=(const Facet& fct);
	bool operator==(const Facet& fct) const;

	//Facet(const Facet&& fct) = delete;
	//Facet& operator=(const Facet&& fct) = delete;

	void operator()(const point& v_r, point& v_Grv);
	void operator()(const point& v_r, atmic_double& g);

protected:
	bool   _initialized; // is initialized?
	int    _id;    // facet ID
	ptvec  _pts;   // points
	ptvec  _L;     // array of vector l vectors; l = pts[i+1] - pts[i]	
	ptvec  _mi;    // array of unit vector mi vectors; mi = Unit(pts[i+1] - pts[i])	
	ptvec  _ni;    // arrys of unit ni vectors; ni = _mi[i] / _n;
	point  _n;     // facet unit normal vector
	valvec _len;   // arry of side lengths
	//valvec _g;     // results gx, gy, gz, gxx, gyy, gzz, gxy, gxz, gyz

	// physical parameters
	// under bodyCCW is understand when facet is seen CCW from outside of the body
	double  _density;		// density of pBodyCCW
	point   _densGrad;		// density gradien of bodyCCW
	bool    _bLin;			// if compute withLinear Density

	// bodyCW -> oposit body, facet is ordered CW
	double  _densityOpos;	// density of oposit body
	point   _densGradOpos;	// density gradien of bodyCCW
	bool    _bLinOpos;		// if compute with linear density

	// facet sign for updating
	double  _sign;		// m_dSign == -1 to remove field; m_dSign == 1 to add field

	// temporary variables
	point   _tmpFld;			// field from body which facet is seen from outside CCW
	point   _tmpFldGrd;		// gradients
	point   _tmpFldOpos;		// oposit facet


public:
	void Init();
	void Init(ptvec& pts);
	void Init(ptvec& pts, double densityCCW, double densityCW = 0.0);

	//double GetSign()          { return _sign; }
	//void SetSign(double sign) { _sign = sign; }

	void Fld_G(const point &v_r, point& v_Grv);
	void Fld_Gz(const point &v_r, atmic_double& g);

	ptvec& Data() { return _pts;}


protected:
	void FldVlado(const point &v_r, double& f);

};

using facet_vec = std::vector < Facet >;


}; // namespace pfld

#endif  // FACET_H_