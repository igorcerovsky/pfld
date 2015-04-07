#include "facet.hpp"

namespace pfld {

Facet::Facet() : _n(0,0,0),
_density(1000),
_densityOpos(1000),
_densGrad(0,0,0),
_densGradOpos(0,0,0),
_sign(1.0),
_mi(),
_ni(),
_pts(),
_L(),
_len()
{
}

Facet::Facet(const Facet& fct) : _id(fct._id),
_n(fct._n),
_density(fct._density),
_densityOpos(fct._densityOpos),
_densGrad(fct._densGrad),
_densGradOpos(fct._densGradOpos),
_sign(fct._sign),
_mi(fct._mi),
_ni(fct._ni)
{
	_pts.assign(fct._pts.begin(), fct._pts.end());
	_L.assign(fct._L.begin(), fct._L.end());
	_len.assign(fct._len.begin(), fct._len.end());
}

Facet& Facet::operator=(const Facet& fct)
{
	_id = fct._id;
	_n = fct._n;
	_density = fct._density;
	_densityOpos = fct._densityOpos;
	_densGrad = fct._densGrad;
	_densGradOpos = fct._densGradOpos;
	_sign = fct._sign;
	_pts.assign(fct._pts.begin(), fct._pts.end());
	_L.assign(fct._L.begin(), fct._L.end());
	_mi = fct._mi;
	_ni = fct._ni;
	_n = fct._n;
	_len.assign(fct._len.begin(), fct._len.end());
	return *this;
}

bool Facet::operator==(const Facet& fct) const
{
	return _n == fct._n;
}

void Facet::operator()(point& v_r, point& v_Grv)
{
	if (!_initialized)
		Init();
	
	FldVlado(v_r, v_Grv);
}

void Facet::Init()
{
	//TRACE("\n Facet Init.\n");
	// side normal vector n; only if the angle between l_th and (l-1)th edge is smaller then PI (180 deg), what is true for triangle 
	//v_n = ( pts[2] - pts[1] ) / ( pts[1] - pts[0] );	// side normal vector (operator /  is vector multiplication see class CPoint 3D)
	_n = (_pts[0] - _pts[1]) / (_pts[1] - _pts[2]);
	_n.Unit();
	//// TRACE("Side normal vector is (%f, %f, %f)\n", (float) v_n.x, (float) v_n.y, (float) v_n.z);

	// for each edge
	const int n(_pts.size());
	_mi.resize(n);
	_ni.resize(n);
	_L.resize(n);
	_len.resize(n);
	for (int i = 0; i < n; i++) {
		//TRACE("   x=%6.1f, y=%6.1f z=%6.1f \n", (float) pts[i].x, (float) pts[i].y, (float) pts[i].z);
		// (a) length and unit vector v_mi	
		if (i != n - 1) {
			_mi[i] = _pts[i + 1] - _pts[i];
			_L[i] = _mi[i];
			_len[i] = _mi[i].Abs();				// edge length
			_mi[i].Unit();
		}
		else {
			_mi[i] = _pts[0] - _pts[i];
			_L[i] = _mi[i];
			_len[i] = _mi[i].Abs();				// edge length
			_mi[i].Unit();
		}
		// TRACE("   length is   %f\n", (float) len[i]);
		// TRACE("   v_mi vector is (%f, %f, %f)\n", (float) v_mi[i].x, (float) v_mi[i].y, (float) v_mi[i].z);

		// edge unit vector v_ni
		_ni[i] = _mi[i] / _n;
		//// TRACE("   v_ni vector is (%f, %f, %f)\n", (float) pt.x, (float) pt.y, (float) pt.z);
	}

	// what kind of computation?
	_bLin = _densGrad.IsZero();
	_bLinOpos = _densGradOpos.IsZero();

	_initialized = true;
}

void Facet::Init(ptvec& pts)
{
	_pts.assign(pts.begin(), pts.end());
	Init();
}

void Facet::Init(ptvec& pts, double densityCCW /*= 1000.0*/, double densityCW /*= 0.0*/)
{
	_pts.assign(pts.begin(), pts.end());
	_density = densityCCW;
	_densityOpos = densityCW;
	Init();
}

// !!! Init MUST be called first
// v_r point distance vector
// v_Grv otput field 
void Facet::FldVlado(point &v_r, point &v_Grv)
{
	const int n( _pts.size() );
	// second part of algorithm, computing field 
	double z, u, v, w, W2, W, U, V, T, f = 0.0, ff = 0.0;

	z = fabs(_n * (_pts.at(0) - v_r));

	for (int i = 0; i < n; i++) {
		// TRACE("p1.x = %6.4f, p1.y = %6.4f, p1.z = %6.4f\n", (float) (pts[i]).x, (float) (pts[i]).y, (float) (pts[i]).z);
		u = _mi[i] * (_pts[i] - v_r);
		v = u + _len[i];
		w = _ni[i] * (_pts[i] - v_r);

		//z = z + EPS;
		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		W = sqrt(W2);
		T = U + V;
		f += w * (sign(v)*log((V + fabs(v)) / W) - sign(u)*log((U + fabs(u)) / W)) -
			2 * z*atan((2 * w*_len[i]) / ((T + _len[i])*fabs(T - _len[i]) + 2 * T*z));
	}
	f *= GRAVCONST;

	v_Grv += _n*f;
	// TRACE("Vlado returns f = %6.12f\n", (float) f * v_n.z);
}


}; // namespace pfld
