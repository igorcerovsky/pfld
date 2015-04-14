#include "facet.hpp"

namespace pfld {

Facet::Facet() : _id(0),
	_sz(0),
	_initialized(false),
	_n(0, 0, 0),
	_mi(),
	_ni(),
	_pts(),
	_L(),
	_len()
{
}

Facet::Facet(const ptvec& pts) : Facet()
{
	_pts.assign(pts.begin(), pts.end());
	_initialized = false;
}

Facet::Facet(const Facet& fct) : _id(fct._id),
	_sz(fct._sz),
	_initialized(fct._initialized),
	_n(fct._n),
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
	_initialized = fct._initialized;
	_sz = fct._sz;
	_n = fct._n;
	_pts.assign(fct._pts.begin(), fct._pts.end());
	_L.assign(fct._L.begin(), fct._L.end());
	_mi.assign(fct._mi.begin(), fct._mi.end());
	_ni.assign(fct._ni.begin(), fct._ni.end());
	_len.assign(fct._len.begin(), fct._len.end());
	return *this;
}

bool Facet::operator==(const Facet& fct) const
{
	return _n == fct._n;
}

void Facet::Init()
{
	if (_initialized)
		return;
	// side normal vector n; only if the angle between l_th and (l-1)th edge is smaller then PI (180 deg), what is true for triangle 
	point::Cross(_pts[0] - _pts[1], _pts[1] - _pts[2], _n);
	_n.Unit();

	// for each edge
	_sz = _pts.size();
	_mi.resize(_sz);
	_ni.resize(_sz);
	_L.resize(_sz);
	_len.resize(_sz);
	int i = 0;
	for (; i < _sz - 1; ++i)
	{
		point::Sub(_pts[i + 1], _pts[i], _mi[i]);
		_L[i] = _mi[i];
		_len[i] = _mi[i].Abs(); // edge length
		_mi[i].Unit();
		 point::Cross(_mi[i], _n, _ni[i]);
	}
	point::Sub(_pts[0], _pts[i], _mi[i]);
	_L[i] = _mi[i];
	_len[i] = _mi[i].Abs();
	_mi[i].Unit();
	point::Cross(_mi[i], _n, _ni[i]);

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
	Init();
}

// !!! Init MUST be called first
// v_r point distance vector
// v_Grv otput field 
void Facet::FldVlado(const point& v_r, double& f)
{
	// second part of the algorithm, computing field 
	double z, u, v, w, W2, W, U, V, T;

	z = fabs(_n * (_pts.at(0) - v_r)) + EPS;

	for (unsigned i = 0; i < _sz; i++)
	{
		point tmp;
		point::Sub(_pts[i], v_r, tmp);
		u = _mi[i] * tmp;
		w = _ni[i] * tmp;
		v = u + _len[i];

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		W = sqrt(W2);
		T = U + V;
		f += w * (sign(v)*log((V + fabs(v)) / W) - sign(u)*log((U + fabs(u)) / W)) -
			2 * z*atan((2 * w*_len[i]) / ((T + _len[i])*fabs(T - _len[i]) + 2 * T*z));
	}
	f *= GRAVCONST;
}


// gravity field for variable density
// ro density gradient vector
// ro0 density in origin of coordinate system
void Facet::Fld_G(const point& v_r, point& v_Grv, const point& ro, const double& ro0)
{
	int n = 3;
	// second part of algorithm, computing field 
	double z, u, v, w, W2, W, U, V, T, L, A, Fi, Fi2, d, Z;
	point f;
	const double ro_r{ ro0 + ro*v_r };

	Z = _n * (_pts[0] - v_r);
	z = fabs(Z) + EPS;

	const double ronz = ro*_n*Z;
	for (int i = 0; i < n; i++) 
	{
		point ptTmp1;
		point::Sub(_pts[i],v_r, ptTmp1);
		u = _mi[i] * ptTmp1;
		v = u + _len[i];
		w = _ni[i] * ptTmp1;

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		W = sqrt(W2);
		T = U + V;
		d = _len[i];
		A = -atan((2 * w*d) / ((T + d)*fabs(T - d) + 2 * T*z));
		if (sign(u) == sign(v)) {
			L = sign(v)*log((V + fabs(v)) / (U + fabs(u)));
		}
		else {
			L = log((V + fabs(v))*(U + fabs(u)) / (W*W));
		}
		Fi = w*L + 2 * z*A;
		Fi2 = d / 4 * ((v + u)*(v + u) / T + T) + W*W*L / 2;
		f += _n*(Fi*(ro_r + ronz) + ro*_ni[i] * Fi2) - ro*(Fi*Z / 2);
	}
	f = f*GRAVCONST;
	v_Grv += f;
}

void Facet::Fld_G(const point& v_r, point& v_Grv)
{
	double f = 0.0;
	FldVlado(v_r, f);
	v_Grv += _n*f;
}


void Facet::Fld_Gz(const point& v_r, double_pfld& g)
{
	double f = 0.0;
	FldVlado(v_r, f);
	f *= _n.z;
	g = g + f;
}

void Facet::operator()(const point& v_r, point& v_Grv)
{
	Fld_G(v_r, v_Grv);
}

void Facet::operator()(const point& v_r, double_pfld& g)
{
	Fld_Gz(v_r, g);
}


}; // namespace pfld
