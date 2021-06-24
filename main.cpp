#include <iostream>
#include <tms-nets.hpp>
#include <chrono>


#include "tms-nets/analysis/analysis.hpp"




void run_sobol_t_tests_for_nets_with_smallest_t_estimate(int s_start = 3, int s_step = 1, int s_end = 15, int m_start = 5, int m_step = 1, int m_end = 25)
{
	auto start = std::chrono::high_resolution_clock::now();
	auto stop  = std::chrono::high_resolution_clock::now();
	
	for (int s = s_start; s < s_end; ++s)
	{
		std::cout << "\n\ns = " << s << ":\n";
		for (int m = m_start; m < m_end; ++m)
		{
			std::cout << "m = " << m << ":\n";
			
			tms::Sobol net(m, s);
			
			std::cout << "POLY: ";
			start = std::chrono::high_resolution_clock::now();
			std::cout << tms::analysis::sobol_t(net, tms::analysis::SobolTestType::poly) << " ";
			stop  = std::chrono::high_resolution_clock::now();
			std::cout << "(" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms)\n";
			
			std::cout << "LIP:  ";
			start = std::chrono::high_resolution_clock::now();
			std::cout << tms::analysis::sobol_t(net, tms::analysis::SobolTestType::lip) << " ";
			stop  = std::chrono::high_resolution_clock::now();
			std::cout << "(" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms)\n";
			
			std::cout << "LIPM: ";
			start = std::chrono::high_resolution_clock::now();
			std::cout << tms::analysis::sobol_t(net, tms::analysis::SobolTestType::lip_mod) << " ";
			stop  = std::chrono::high_resolution_clock::now();
			std::cout << "(" << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms)\n";
		}
	}
}



int main(void)
{
	run_sobol_t_tests_for_nets_with_smallest_t_estimate();

	return 0;
}