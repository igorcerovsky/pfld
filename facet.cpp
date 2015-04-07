#include "facet.hpp"

namespace pfld {
	Facet::Facet()
	{

	}

	Facet::Facet(const Facet& fct)
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
		return true;
	}


	void Facet::Init()
	{

	}

	void Facet::Init(fctpointvector& pts)
	{
	}

	void Facet::Init(fctpointvector& pts, double densityCCW /*= 1000.0*/, double densityCW /*= 0.0*/)
	{

	}

	void Facet::FldVlado(fctpoint &v_r, fctpoint &v_Grv)
	{

	}


}; // namespace pfld
