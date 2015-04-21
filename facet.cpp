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
	size_t i = 0;
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

///////////////////////////////////////////////////////////////////////////////
//********** field Vlado Pohanka **********************************************

// !!! Init MUST be called first
// r point distance vector
// grv otput field 
void Facet::FldVlado(const point& r, double& f)
{
	// second part of the algorithm, computing field 
	double z, u, v, w, W2, W, U, V, T;

	z = fabs(_n * (_pts.at(0) - r)) + EPS;

	for (unsigned i = 0; i < _sz; i++)
	{
		point tmp;
		point::Sub(_pts[i], r, tmp);
		u = _mi[i] * tmp;
		w = _ni[i] * tmp;
		v = u + _len[i];

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		W = sqrt(W2);
		T = U + V;
		f += w * (sign(v)*log((V + fabs(v)) / W) - sign(u)*log((U + fabs(u)) / W)) -
			2.0 * z*atan((2.0 * w *_len[i]) / ((T + _len[i])*fabs(T - _len[i]) + 2.0 * T*z));
	}
	f *= GRAVCONST;
}


// gravity field for linearly variable density
// ro density gradient vector
// ro0 density in origin of coordinate system
void Facet::Fld_Gz(const point& r, const point& ro, const double& ro0, double& gz)
{
	// second part of algorithm, computing field 
	double f{ 0.0 };
	const double ro_r{ ro0 + ro*r };

	double Z{ _n * (_pts[0] - r) };
	double z{ fabs(Z) + EPS };

	const double ronz{ ro*_n*Z };
	for (size_t i = 0; i < _sz; ++i)
	{
		double u, v, w, W2, U, V, T, L, A, Fi, Fi2;
		point ptTmp1;
		point::Sub(_pts[i], r, ptTmp1);
		u = _mi[i] * ptTmp1;
		v = u + _len[i];
		w = _ni[i] * ptTmp1;

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		T = U + V;
		double& d = _len[i];
		A = -atan((2.0 * w*d) / ((T + d)*fabs(T - d) + 2.0 * T*z));
		if (sign(u) == sign(v)) {
			L = sign(v)*log((V + fabs(v)) / (U + fabs(u)));
		}
		else {
			L = log((V + fabs(v))*(U + fabs(u)) / W2);
		}
		Fi = w*L + 2.0 * z*A;
		Fi2 = d * 0.25 * ((v + u)*(v + u) / T + T) + W2*L * 0.5;
		f += _n.z*(Fi*(ro_r + ronz) + ro*_ni[i] * Fi2) - ro.z*(Fi*Z * 0.5);
	}
	f = f*GRAVCONST;
	gz += f;
}


// gravity field for linearly variable density
// ro density gradient vector
// ro0 density in origin of coordinate system
void Facet::Fld_G(const point& r, const point& ro, const double& ro0, point& grv)
{
	// second part of algorithm, computing field 
	point f;
	const double ro_r{ ro0 + ro*r };

	double Z{ _n * (_pts[0] - r) };
	double z{ fabs(Z) + EPS };

	const double ronz = ro*_n*Z;
	for (size_t i = 0; i < _sz; ++i) 
	{
		double u, v, w, W2, U, V, T, L, A, Fi, Fi2;
		point ptTmp1;
		point::Sub(_pts[i],r, ptTmp1);
		u = _mi[i] * ptTmp1;
		v = u + _len[i];
		w = _ni[i] * ptTmp1;

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		T = U + V;
		double& d = _len[i];
		A = -atan((2.0 * w*d) / ((T + d)*fabs(T - d) + 2.0 * T*z));
		if (sign(u) == sign(v)) {
			L = sign(v)*log((V + fabs(v)) / (U + fabs(u)));
		}
		else {
			L = log((V + fabs(v))*(U + fabs(u)) / W2);
		}
		Fi = w*L + 2.0 * z*A;
		Fi2 = d * 0.25 * ((v + u)*(v + u) / T + T) + W2*L * 0.5;
		f += _n*(Fi*(ro_r + ronz) + ro*_ni[i] * Fi2) - ro*(Fi*Z *0.5);
	}
	f = f*GRAVCONST;
	grv += f;
}

void Facet::Fld_G(const point& r, point& grv)
{
	double f = 0.0;
	FldVlado(r, f);
	grv += _n*f;
}


void Facet::Fld_Gz(const point& r, double_pfld& g)
{
	double f = 0.0;
	FldVlado(r, f);
	f *= _n.z;
	g = g + f;
}

void Facet::operator()(const point& r, point& grv)
{
	Fld_G(r, grv);
}

void Facet::operator()(const point& r, double_pfld& g)
{
	Fld_Gz(r, g);
}

//********** end field Vlado Pohanka ******************************************
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//********** field Singh && Guptasarma ****************************************
// measuring point MUST be coordinate system origin => first transform points
// v_M			magnetization vector
// v_Mag		output magnetic field
// r			measuring point
void Facet::FldGS(const point& r, const point v_M, point& mag, point& grv)
{
	// shifted points, make local copy of shifted points, do not change origanal points
	ptvec	spts(_sz);	
	point	shf{ r * (-1) };
	for (size_t i = 0; i<_sz; i++) {
		point::Add(_pts[i], shf, spts[i]);
	}

	double P = 0., Q = 0., R = 0.;
	double dOmega = SolidAngle(spts);

	// for all edges compute length, P, Q, R
	for (size_t i = 0; i<_sz; i++) {
		double r = spts[i].Abs();
		double b = 2 * (spts[i] * _L[i]);
		double h = r + b / (2 * _len[i]);
		double I;
		if (h > EPS)
			I = (1.0 / _len[i])*log((sqrt(_len[i]*_len[i] + b + r*r) + _len[i] + b / (2.0 * _len[i])) / h);
		else
			I = (1.0 / _len[i]) * log(fabs(_len[i] - r) / r);

		P += I*_L[i].x;
		Q += I*_L[i].y;
		R += I*_L[i].z;
	}

	point f{ dOmega*_n.x + Q*_n.z - R*_n.y,
		dOmega*_n.y + R*_n.x - P*_n.z,
		dOmega*_n.z + P*_n.y - Q*_n.x };

	// magnetic field
	double s = v_M * _n; // surface pole density
	mag += f * s;

	// gravity field
	double d = spts[0] * _n;
	grv += f * d * GRAVCONST;	// kappa
}


// computes solid angle as seen from the ORIGIN of coordinate system
// pts   array  of pointers to polygon points
double Facet::SolidAngle(ptvec& pts)
{
	double Omega = 0.;	// solid angle
	// is the polygon seen from the otside or inside?
	double dFi = 0.;
	double dInOut{ _n * pts[1] };
	if (dInOut == 0.)
		return 0.;

	std::vector<point*> pp(_sz + 2);
	auto itPP = pp.begin();
	*itPP = &(*(--pts.end()));
	itPP++;
	for (auto it = pts.begin(); it != pts.end(); ++it, ++itPP)
		*itPP = &(*it);
	*itPP = &(*(pts.begin()));
	if (dInOut >  0.) {
		std::reverse(pp.begin(), pp.end());
	}

	for (size_t i = 0; i<_sz; i++) {
		point n1, n2;
		point::Cross(*pp[i + 1], *pp[i], n1);
		n1.Unit();
		point::Cross(*pp[i + 2], *pp[i + 1], n2);
		n2.Unit();
		double dPerp = *pp[i+2] * n1;
		double b = n1 * n2;
		if (b<-1.0) { b = -1.0; }
		if (b>1.0)  { b = 1.0; }
		double a = PI - acos(b);
		if (dPerp < 0.) {
			a = PI2 - a;
		}
		dFi += a;
	}
	Omega = dFi - (_sz - 2.0)*PI;
	if (dInOut > 0)
		Omega = -Omega;

	return Omega;
}
//********** end field Singh && Guptasarma ************************************
///////////////////////////////////////////////////////////////////////////////




}; // namespace pfld
