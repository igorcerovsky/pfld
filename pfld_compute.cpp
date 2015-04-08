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
	unsigned long const length = std::distance(first, last);

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


void Field_Gz(pfld::facet_vec& facets, pfld::ptvec& field_points, pfld::valvec& out_field)
{
	parallel_init(facets.begin(), facets.end());
}

}; // namespace pfld