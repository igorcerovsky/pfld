#include "pfld_compute.hpp"
#include <thread>


namespace pfld {

	
void Field_G(pfld::facet_vec& facets, pfld::ptvec& field_points, pfld::ptvec& out_field)
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


template<typename Iterator>
struct init_block
{
	void operator()(Iterator first, Iterator last)
	{
		auto it = first;
		for (; it != last; ++it)
			it->Init();
	}
};

template<typename Iterator>
void parallel_init(Iterator first, Iterator last)
{
	unsigned long const length = unsigned long const(std::distance(first, last));

	if (length==0)
		return;

	unsigned long const min_per_thread = 25;
	unsigned long const max_threads =
		(length + min_per_thread - 1) / min_per_thread;

	unsigned long const hardware_threads = std::thread::hardware_concurrency();

	auto const num_threads =
		std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);

	auto const block_size = length / num_threads;

	std::vector<std::thread> threads(num_threads);

	Iterator block_start = first;
	for (unsigned i = 0; i<num_threads; ++i)
	{
		Iterator block_end = block_start;
		std::advance(block_end, block_size);
		threads[i] = std::thread(
			init_block<Iterator>(), block_start, block_end);
		block_start = block_end;
	}

	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}


template<typename IteratorPoints, typename Facets, typename IteratorField>
struct compute_block
{
	void operator()(IteratorPoints first, IteratorPoints last, Facets& facets, IteratorField itFld)
	{
		for (auto it = first; it != last; ++it, ++itFld)
		{
			point& pt = *it;
			double_pfld f(0.0);
			for (auto itFcs = facets.begin(); itFcs != facets.end(); ++itFcs ) {
				itFcs->Fld_Gz(pt, f);
			}
			*itFld = f;
		}
	}
};

template<typename IteratorPoints, typename Facets, typename IteratorField>
void parallel_compute(IteratorPoints first, IteratorPoints last, Facets& facets, IteratorField itFld)
{
	unsigned long const length = unsigned long const(std::distance(first, last));

	if (length == 0)
		return;

	unsigned long const min_per_thread = 25;
	unsigned long const max_threads =
		(length + min_per_thread - 1) / min_per_thread;

	unsigned long const hardware_threads = std::thread::hardware_concurrency();

	auto const num_threads =
		std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);

	auto const block_size = length / num_threads;

	std::vector<std::thread> threads(num_threads);

	IteratorPoints block_start = first;
	for (unsigned i = 0; i<num_threads; ++i)
	{
		auto block_end = block_start;
		std::advance(block_end, block_size);
		threads[i] = std::thread(
			compute_block<IteratorPoints, Facets, IteratorField>(),
			block_start, block_end, facets, itFld);
		block_start = block_end;
		std::advance(itFld, block_size);
	}

	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}


// parallell approach
void Field_Gz(pfld::facet_vec& facets, pfld::ptvec& fldPoints, pfld::valvec& outFld)
{
	parallel_init(facets.begin(), facets.end());

	parallel_compute(fldPoints.begin(), fldPoints.end(), facets, outFld.begin());
}


// naive approach
void Field_Gz__(pfld::facet_vec& facets, pfld::ptvec& fldPoints, pfld::valvec& outFld)
{
	for (auto it = facets.begin(); it != facets.end(); ++it)
		it->Init();

	auto itFld = outFld.begin();
	for (auto it = fldPoints.begin(); it != fldPoints.end(); ++it, ++itFld)
	{
		point& pt = *it;
		double_pfld f(0.0);
		for (auto itFcs = facets.begin(); itFcs != facets.end(); ++itFcs)
		{
			itFcs->Fld_Gz(pt, f);
		}
		*itFld = f;
	}

}

}; // namespace pfld