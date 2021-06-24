#ifndef TMS_NETS_ANALYSIS_SOBOL_HPP
#define TMS_NETS_ANALYSIS_SOBOL_HPP

#include "../sobol.hpp"

namespace tms::analysis
{
	
	enum SobolTestType
	{
		poly,
		lip,
		lip_mod
	};
	
	/** Computes t-value of Sobol' net over GF[2] with the certain method
	 * @param net – Sobol' net
	 * @param test_type – id of computation method
	 *
	 */
	BasicInt find_sobol_defect(Sobol const &net, SobolTestType test_type, int t_estimate = -1);
	
};


#endif
