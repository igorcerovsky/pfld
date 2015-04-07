// pfld.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "pfld.hpp"
#include "facet.hpp"

PFLD_API void Field_G(double* facets, const int num_facets,
	double* field_points, double* out_field, const int num_points)
{
	// TODO
}

PFLD_API void Field_G(pfld::facet_vec& facets, pfld::ptvec& field_points, pfld::ptvec& out_field)
{
	if (field_points.size() != out_field.size())
		return;
	pfld::point pt(0.0, 0.0, 0.0);
	pfld::point g(0.0, 0.0, 0.0);
	auto itf = out_field.begin();
	for (auto itp = field_points.begin(); itp != field_points.end(); ++itp, ++itf)
	{
		for (auto it = facets.begin(); it != facets.end(); ++it)
			it->operator()(*itp, *itf);
	}
}