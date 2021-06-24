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
	
	/** Returns t-value of Sobol' net over GF[2] evaluated with the certain algorithm
	 * @param net – Sobol' net
	 * @param test_type – enum value of applied algorithm method
	 * @param t_estimate – (optional) t-value estimate for the net
	 */
	BasicInt sobol_t(Sobol const &net, SobolTestType test_type, int t_estimate = -1);
	
};


#endif
