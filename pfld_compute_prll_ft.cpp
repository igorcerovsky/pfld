#include "pfld_compute.hpp"
#include <thread>
#include <future>
#include <utility> 

namespace pfld {


class join_threads
{
	std::vector<std::thread>& threads;
public:
	explicit join_threads(std::vector<std::thread>& threads_) :
		threads(threads_)
	{}
	~join_threads()
	{
		for (unsigned long i = 0; i < threads.size(); ++i)
		{
			if (threads[i].joinable())
				threads[i].join();
		}
	}
};

template<typename Iter>
struct init_block_ft
{
	int operator()(Iter first, Iter last)
	{
		for (auto it = first; it != last; ++it)
		{
			it->Init();
		}
		return int(std::distance(first, last));
	}
};

template<typename Iter>
void parallel_init_ft(Iter first, Iter last)
{
	unsigned long const length = (unsigned long const)std::distance(first, last);

	if (!length)
		return;

	unsigned long const min_per_thread = 25;
	unsigned long const max_threads =
		(length + min_per_thread - 1) / min_per_thread;

	unsigned long const hardware_threads =
		std::thread::hardware_concurrency();

	unsigned long const num_threads =
		std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);

	const auto actual_threads = num_threads - 1;
	unsigned long const block_size = length / actual_threads;

	std::vector<std::future<int>> futures(actual_threads);
	std::vector<std::thread> threads(actual_threads);
	join_threads joiner(threads);

	Iter block_start = first;
	Iter block_end = block_start;
	for (unsigned long i = 0; i < actual_threads; ++i)
	{
		block_end = block_start;
		if (i == (actual_threads - 1))
			block_end = last;
		else
			std::advance(block_end, block_size);

		std::packaged_task<int(Iter, Iter)> task((init_block_ft<Iter>()));
		futures[i] = task.get_future();
		threads[i] = std::thread(std::move(task), block_start, block_end);
		block_start = block_end;
	}
}


template<typename ItPts, typename Facets, typename FldType>
struct cmp_block_ft
{
	FldType operator()(ItPts first, ItPts last, Facets& facets)
	{
		FldType fld(std::distance(first, last), 0.);
		auto itFld = fld.begin();
		for (auto it = first; it != last; ++it, ++itFld)
		{
			double_pfld f(0.0);
			for (auto itFcs = facets.begin(); itFcs != facets.end(); ++itFcs) {
				itFcs->Fld_Gz(*it, f);
			}
			*itFld = f;
		}
		return fld;
	}
};

template<typename ItPts, typename Facets, typename OutField>
void parallel_cmp_ft(ItPts first, ItPts last, Facets& facets, OutField& field)
{
	unsigned long const length = (unsigned long const)std::distance(first, last);

	if (!length)
		return;

	unsigned long const min_per_thread = 25;
	unsigned long const max_threads =
		(length + min_per_thread - 1) / min_per_thread;

	unsigned long const hardware_threads =
		std::thread::hardware_concurrency();

	unsigned long const num_threads =
		std::min(hardware_threads != 0 ? hardware_threads : 2, max_threads);

	const auto actual_threads = num_threads;
	unsigned long const block_size = length / actual_threads;

	std::vector<std::future<valvec>> futures(actual_threads);
	std::vector<std::thread> threads(actual_threads);
	join_threads joiner(threads);

	ItPts block_start = first;
	ItPts block_end = block_start;
	for (unsigned long i = 0; i < actual_threads; ++i)
	{
		block_end = block_start;
		if (i == (actual_threads - 1))
			block_end = last;
		else
			std::advance(block_end, block_size);

		std::packaged_task<valvec(ItPts, ItPts, Facets&)>
			task((cmp_block_ft<ItPts, Facets, valvec>()));
		futures[i] = task.get_future();
		threads[i] = std::thread(std::move(task), block_start, block_end, facets);
		block_start = block_end;
	}
	for (unsigned long i = 0; i < actual_threads; ++i)
	{
		valvec thread_result = std::move(futures[i].get());
		field.insert(field.end(), thread_result.begin(), thread_result.end());
	}
}

// parallell approach with futute
// output field is added to the outFld vector in order of fldPoints
void Field_Gz_(pfld::facet_vec& facets, pfld::ptvec& fldPoints, pfld::valvec& outFld)
{
	parallel_init_ft(facets.begin(), facets.end());

	parallel_cmp_ft(fldPoints.begin(), fldPoints.end(), facets, outFld);
}


} //namespace pfld
