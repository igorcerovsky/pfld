#include "facet_t.hpp"

namespace pfld {

template<typename T>
Facet<T>::Facet(void) : _id(-1),
	_sz(0),
	_initialized(false),
	_n(),
	_mi(),
	_ni(),
	_pts(),
	_L(),
	_len()
{
}

template<typename T>
Facet<T>::Facet(const ptvec& pts) : Facet()
{
	_pts.assign(pts.begin(), pts.end());
	_initialized = false;
}

template<typename T>
Facet<T>::Facet(const Facet& fct) : _id(fct._id),
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

template<typename T>
Facet<T>& Facet<T>::operator=(const Facet<T>& fct)
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

template<typename T>
bool Facet<T>::operator==(const Facet& fct) const
{
	return _n == fct._n;
}

template<typename T>
void Facet<T>::Init()
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

template<typename T>
void Facet<T>::Init(ptvec& pts)
{
	_pts.assign(pts.begin(), pts.end());
	Init();
}


// !!! Init MUST be called first
// r point distance vector
// grv otput field 
template<typename T>
void Facet<T>::FldVlado(const point& r, T& f)
{
	// second part of the algorithm, computing field 
	T z, u, v, w, W2, W, U, V, TT;

	z = fabs(_n * (_pts.at(0) - r)) + T(EPS);

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
		TT = U + V;
		f += w * (sign(v)*log((V + fabs(v)) / W) - sign(u)*log((U + fabs(u)) / W)) -
			2 * z*atan((2 * w *_len[i]) / ((TT + _len[i])*fabs(TT - _len[i]) + 2 * TT*z));
	}
	f *= T(GRAVCONST);
}


// gravity field for linearly variable density
// ro density gradient vector
// ro0 density in origin of coordinate system
template<typename T>
void Facet<T>::Fld_Gz(const point& r, const point& ro, const T& ro0, T& gz)
{
	// second part of algorithm, computing field 
	T f{ 0.0 };
	const T ro_r{ ro0 + ro*r };

	T Z{ _n * (_pts[0] - r) };
	T z{ fabs(Z) + EPS };

	const T ronz{ ro*_n*Z };
	for (size_t i = 0; i < _sz; ++i)
	{
		T u, v, w, W2, U, V, T, L, A, Fi, Fi2;
		point ptTmp1;
		point::Sub(_pts[i], r, ptTmp1);
		u = _mi[i] * ptTmp1;
		v = u + _len[i];
		w = _ni[i] * ptTmp1;

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		T = U + V;
		T& d = _len[i];
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
template<typename T>
void Facet<T>::Fld_G(const point& r, const point& ro, const T& ro0, point& grv)
{
	// second part of algorithm, computing field 
	point f;
	const T ro_r{ ro0 + ro*r };

	T Z{ _n * (_pts[0] - r) };
	T z{ fabs(Z) + EPS };

	const T ronz = ro*_n*Z;
	for (size_t i = 0; i < _sz; ++i) 
	{
		T u, v, w, W2, U, V, T, L, A, Fi, Fi2;
		point ptTmp1;
		point::Sub(_pts[i],r, ptTmp1);
		u = _mi[i] * ptTmp1;
		v = u + _len[i];
		w = _ni[i] * ptTmp1;

		W2 = w*w + z*z;
		U = sqrt(u*u + W2);
		V = sqrt(v*v + W2);
		T = U + V;
		T& d = _len[i];
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

template<typename T>
void Facet<T>::Fld_G(const point& r, point& grv)
{
	T f{ 0.0 };
	FldVlado(r, f);
	grv += _n*f;
}

template<typename T>
void Facet<T>::Fld_Gz(const point& r, double_pfld& g)
{
	T f = 0.0;
	FldVlado(r, f);
	f *= _n.z;
	g = g + f;
}

template<typename T>
void Facet<T>::operator()(const point& r, point& grv)
{
	Fld_G(r, grv);
}

template<typename T>
void Facet<T>::operator()(const point& r, double_pfld& g)
{
	Fld_Gz(r, g);
}


// M			magnetization vector
// Mag		output magnetic field
// r			measuring point - oposit direction 
// mag			computed magnetic field
// grv			computed gravity field
template<typename T>
void Facet<T>::FldGS(const point& r, const point M, point& mag, point& grv)
{
	point f;
	FldGS(r, f);

	// magnetic field
	T s = M * _n; // surface pole density
	mag += f * s;

	// gravity field
	T d = (_pts[0] + (r*-1.)) * _n;
	grv += f * d * T(GRAVCONST);	// kappa
}

template<typename T>
void Facet<T>::FldGS_M(const point& r, const point M, point& mag)
{
	point f;
	FldGS(r, f);
	mag += f * (M * _n);
}

template<typename T>
void Facet<T>::FldGS_G(const point& r, point& grv)
{
	point f;
	FldGS(r, f);
	grv += f * ((_pts[0] + (r*-1.)) * _n) * GRAVCONST;	// kappa
}

template<typename T>
void Facet<T>::FldGS_Gz(const point& r, T& gz)
{
	point f;
	FldGS(r, f);
	gz += f.z * ((_pts[0] + (r*-1.)) * _n) * T(GRAVCONST);	// kappa
}

// r  measuring point
// f  output fielkd
template<typename T>
void Facet<T>::FldGS(const point& r, point& f)
{
	// shifted points, make local copy of shifted points, do not change origanal points
	ptvec	spts(_sz);
	point	shf{ r * (-1) };
	for (size_t i = 0; i<_sz; i++) {
		point::Add(_pts[i], shf, spts[i]);
	}

	T P{ 0. }, Q{ 0. }, R{ 0. };
	T dOmega = SolidAngle(spts, _n * spts[1], _sz);

	// for all edges compute length, P, Q, R
	for (size_t i = 0; i<_sz; i++) {
		T r = spts[i].Abs();
		T b = 2 * (spts[i] * _L[i]);
		T h = r + b / (2 * _len[i]);
		T I;
		if (h > EPS)
			I = (1 / _len[i])*log((sqrt(_len[i] * _len[i] + b + r*r) + _len[i] + b / (2 * _len[i])) / h);
		else
			I = (1 / _len[i]) * log(fabs(_len[i] - r) / r);

		P += I*_L[i].x;
		Q += I*_L[i].y;
		R += I*_L[i].z;
	}

	f = point{ dOmega*_n.x + Q*_n.z - R*_n.y,
		dOmega*_n.y + R*_n.x - P*_n.z,
		dOmega*_n.z + P*_n.y - Q*_n.x };
}

// computes solid angle as seen from the ORIGIN of coordinate system
// pts   array  of pointers to polygon points
template<typename T>
T Facet<T>::SolidAngle(ptvec& pts, const T inOut, const size_t sz)
{
	T Omega{ 0. };	// solid angle
	// is the polygon seen from the otside or inside?
	T dFi = 0.;
	if (inOut == 0.)
		return 0.;

	std::vector<point*> pp(sz + 2);
	auto itPP = pp.begin();
	*itPP = &(*(--pts.end()));
	itPP++;
	for (auto it = pts.begin(); it != pts.end(); ++it, ++itPP)
		*itPP = &(*it);
	*itPP = &(*(pts.begin()));
	if (inOut >  0.) {
		std::reverse(pp.begin(), pp.end());
	}

	for (size_t i = 0; i<sz; i++) {
		point n1, n2;
		point::Cross(*pp[i + 1], *pp[i], n1);
		n1.Unit();
		point::Cross(*pp[i + 2], *pp[i + 1], n2);
		n2.Unit();
		T dPerp = *pp[i + 2] * n1;
		T b = n1 * n2;
		if (b<-1.0) { b = -1.0; }
		if (b>1.0)  { b = 1.0; }
		T a = T(PI) - acos(b);
		if (dPerp < 0.) {
			a = T(PI2) - a;
		}
		dFi += a;
	}
	Omega = dFi - (sz - 2)*T(PI);
	if (inOut > 0)
		Omega = -Omega;

	return Omega;
}



}; // namespace pfld
