#ifndef FACET_H_
#define FACET_H_


#include <vector>
#include <array>


namespace pfld {

template <typename T>
class Point3D
{
using Points = std::array< T, 3 >;

public:
	Point3D() {};
	virtual ~Point3D() {};
	const Point3D& operator=(const Point3D& pt)  { _coord = pt._coord; return *this; };
	const Point3D& operator=(const Point3D&& pt) { _coord = pt._coord; return *this; };

protected:
	Points _coord;
};


enum class FacetType { normal, oposite };

class Facet
{
using fctvector = std::vector < double > ;
using fctpoint = Point3D < double >;
using fctpointvector = std::vector < fctpoint > ;

public:
	virtual ~Facet();
	Facet();
	Facet(const Facet& fct);

	Facet& operator=(const Facet& fct);
	bool operator==(const Facet& fct) const;

	Facet(const Facet&& fct) {};
	Facet& operator=(const Facet&& fct) {};

protected:
	int       _id;    // facet ID
	fctvector _pts;   // points
	fctvector _L;     // array of vector l vectors; l = pts[i+1] - pts[i]	
	fctpoint  _mi;    // array of unit vector mi vectors; mi = Unit(pts[i+1] - pts[i])	
	fctpoint  _ni;    // arrys of unit ni vectors; ni = _mi[i] / _n;
	fctpoint  _n;     // facet unit normal vector
	fctvector _len;   // arry of side lengths
	fctvector _g;     // results gx, gy, gz, gxx, gyy, gzz, gxy, gxz, gyz

	// physical parameters
	// under bodyCCW is understand when facet is seen CCW from outside of the body
	double        _density;		// density of pBodyCCW
	fctpoint _densGrad;		// density gradien of bodyCCW
	bool     _bLin;			// if compute withLinear Density

	// bodyCW -> oposit body, facet is ordered CW
	double        _densityOpos;	// density of oposit body
	fctpoint _densGradOpos;	// density gradien of bodyCCW
	bool     _bLinOpos;		// if compute with linear density

	// facet sign for updating
	double        _sign;		// m_dSign == -1 to remove field; m_dSign == 1 to add field

	// temporary variables
	fctpoint  _tmpFld;			// field from body which facet is seen from outside CCW
	fctpoint  _tmpFldGrd;		// gradients
	fctpoint  _tmpFldOpos;		// oposit facet


public:
	void Init();
	void Init(fctpointvector& pts);
	void Init(fctpointvector& pts, double densityCCW=1000.0, double densityCW = 0.0);

	void FldVlado(fctpoint &v_r, fctpoint &v_Grv);

	double GetSign()          { return _sign; }
	void SetSign(double sign) { _sign = sign; }


protected:
	double SolidAngle(fctpoint *spts);

};

}; // namespace pfld

#endif  // FACET_H_