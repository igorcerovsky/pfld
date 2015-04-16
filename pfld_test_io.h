#ifndef _PFLD_TEST_IO_H__
#define _PFLD_TEST_IO_H__

#include "facet.hpp"
#include <string>


namespace pfld {


void GetFacets(pfld::facet_vec& facets, const std::string sfile, const int n, bool bGenerate=false);
void GetFieldPoints(pfld::ptvec& points, const std::string sfile, const int n, bool bGenerate = false);
void SaveResults(const pfld::valvec& res, const std::string sfile);
void LoadResults(pfld::valvec& res, const std::string sfile, const int n);

}


#endif //_PFLD_TEST_IO_H__