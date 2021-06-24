#include "../../include/tms-nets/analysis/analysis.hpp"

namespace
{
	using Dimensions   = std::vector<int>;
	using Composition  = std::vector<int>;
	using Compositions = std::vector<Composition>;
	
	using Callback     = std::function<bool(Dimensions const &, Compositions const &, Compositions::iterator &)>;
	
	
	
	class ClosurePOLY
	{
		std::vector<tms::Polynomial> polys;
		
		int xi_coeffs;
		int upper_xi_degree;
		
		int s;
		int m;
		
		std::vector<int> e;
		
		tms::Polynomial unit_poly;
		tms::Polynomial xi;
		tms::Polynomial beta_u_star;
		
		uintmax_t beta_u_star_coeffs;
		
		
	public:
		
		ClosurePOLY(tms::Sobol const &net) :
		unit_poly(tms::gf2poly::make_gf2poly({1})),
		xi(unit_poly),
		beta_u_star(unit_poly)
		{
			polys = net.polynomials();
			
			xi_coeffs = 1;
			upper_xi_degree = 0;
			
			s = net.s();
			m = net.m();
			
			e = std::vector<int>(s, 0);
			for (int i = 0; i < s; ++i)
			{
				e[i] = static_cast<int>(polys[i].degree());
			}
			
			beta_u_star_coeffs = 1;
		}
		
		
		
		bool operator()(Dimensions const &dims, Compositions const &d_vecs, Compositions::iterator &iter_d_vecs)
		{
			int dim_length = static_cast<int>(dims.size());
			
			bool allright = true;
			
			std::vector<int> prev_u_star(dim_length), u_star(dim_length), r_star(dim_length);
			std::vector<tms::Polynomial> pi_star_inv(dim_length, unit_poly), exp_pi(dim_length, unit_poly);
			
			int e_star = 0;
			for (int dim_i = 0; dim_i < dim_length; ++dim_i)
			{
				u_star[dim_i] = ((*iter_d_vecs)[dim_i] - 1) / e[dims[dim_i]] + 1;
				e_star += e[dims[dim_i]]*u_star[dim_i];
			}
			
			upper_xi_degree = e_star - m;
			
			if ( e_star > m )
			{
				//initialize exp_pi and pi_star_inv:
				for (int dim_i = 0; dim_i < dim_length; ++dim_i)
				{
					r_star[dim_i] = ((*iter_d_vecs)[dim_i] - 1) % e[dims[dim_i]];
					exp_pi[dim_i] = tms::gf2poly::pow(polys[dims[dim_i]], u_star[dim_i]);
				}
				
				for (int dim_i = 0; dim_i < dim_length; ++dim_i)
				{
					pi_star_inv[dim_i] = unit_poly;
					
					for (int dim_ii = 0; dim_ii < dim_i; ++dim_ii)
					{
						pi_star_inv[dim_i] *= exp_pi[dim_ii];
					}
					for (int dim_ii = dim_i + 1; dim_ii < dim_length; ++dim_ii)
					{
						pi_star_inv[dim_i] *= exp_pi[dim_ii];
					}
					
					pi_star_inv[dim_i] = tms::gf2poly::find_inverse(pi_star_inv[dim_i], polys[dims[dim_i]]);
				}
				
				do
				{
					//process current u_star:
					xi_coeffs = 1;
					while ( allright && xi_coeffs < (1ULL << upper_xi_degree) )
					{
						allright = false;
						
						xi = tms::gf2poly::number_to_poly(xi_coeffs);
						
						for (int dim_i = 0; !allright && dim_i < dim_length; ++dim_i)
						{
							beta_u_star = (xi * pi_star_inv[dim_i]) % polys[dims[dim_i]];
							beta_u_star_coeffs = tms::gf2poly::poly_to_number(beta_u_star);
							
							allright = ( (beta_u_star_coeffs & ((1ULL << (e[dims[dim_i]] - 1 - r_star[dim_i])) - 1)) != 0 );
						}
						
						++xi_coeffs;
					}
					
					++iter_d_vecs;
					
					for (int dim_i = 0; iter_d_vecs != d_vecs.end() && dim_i < dim_length; ++dim_i)
					{
						prev_u_star[dim_i] = u_star[dim_i];
						u_star[dim_i]      = ((*iter_d_vecs)[dim_i] - 1) / e[dims[dim_i]] + 1;
						r_star[dim_i]      = ((*iter_d_vecs)[dim_i] - 1) % e[dims[dim_i]];
					}
				}
				while ( allright && prev_u_star == u_star && iter_d_vecs != d_vecs.end() );
				
			}//if ( !optimized || e_star > m )
			else
			{
				++iter_d_vecs;
			}
			
			return allright;
		}
		
	};
	
	
	class ClosureLIP
	{
		int m;
		int s;
		
		std::vector<int> e;
		
		std::vector< std::vector<uintmax_t> > gen_mats_row_numbers;
		
		bool optimized;
		
		
	public:
		
		ClosureLIP(tms::Sobol const &net, bool flag = false) :
		m(net.m()),
		s(net.s()),
		e(s, 0),
		gen_mats_row_numbers(s, std::vector<uintmax_t>(m, 0)),
		optimized(flag)
		{
			for (int i = 0; i < s; ++i)
			{
				e[i] = static_cast<int>(net.polynomial(i).degree());
				tms::GenNum gen_num = net.generating_numbers(i);
				for (int j = 0; j < m; ++j)
				{
					for (int k = 0; k < m; ++k)
					{
						gen_mats_row_numbers[i][j] |= (static_cast<uintmax_t>(gen_num.get_bit(j, k)) << (m - 1 - k));
					}
				}
			}
		}
		
		
		
		bool operator()(Dimensions const &dims, Compositions const &d_vecs, Compositions::iterator &iter_d_vecs)
		{
			int dim_length = static_cast<int>(dims.size());
			
			bool allright = true;
			
			int e_star = 0;
			for (int dim_i = 0; dim_i < dim_length; ++dim_i)
			{
				e_star += e[dims[dim_i]]*(((*iter_d_vecs)[dim_i] - 1) / e[dims[dim_i]] + 1);
			}
			
			if ( !optimized || e_star > m )
			{
				std::map<int, uintmax_t> composition_summary;
				
				auto log2_int = [](uintmax_t arg) {
					int degree = 0;
					while ( arg >>= 1 )
					{ ++degree; }
					return degree;
				};
				
				std::function<bool(int, uintmax_t)> process_insertion = [&](int key, uintmax_t value) {
					
					bool not_inserted = true, not_zero = value != 0;
					
					while ( not_zero && not_inserted )
					{
						auto iter_summary = composition_summary.find(key);
						
						if ( iter_summary == composition_summary.end() )
						{
							composition_summary.insert(std::make_pair(key, value));
							not_inserted = false;
							not_zero     = value != 0;
						}
						else
						{
							value   ^= iter_summary->second;
							key      = log2_int(value);
							not_zero = value != 0;
						}
					}
					
					return not_zero;
				};
				
				
				for (int dim_i = 0; allright && dim_i < dim_length; ++dim_i)
				{
					for (int d_vec_i = 0; allright && d_vec_i < (*iter_d_vecs)[dim_i]; ++d_vec_i)
					{
						uintmax_t cur_row_number = gen_mats_row_numbers[dims[dim_i]][d_vec_i];
						allright = process_insertion(log2_int(cur_row_number), cur_row_number);
					}
				}
			}
			
			++iter_d_vecs;
			
			return allright;
		}
	};
	
	
	bool
	get_next_composition(Composition &comp, unsigned sum_of_elements)
	{
		int local_sum = 0;
		int i = 0, length = static_cast<int>(comp.size());
		
		while ( i < length - 1 && local_sum < (sum_of_elements - length + i) )
		{
			local_sum += comp[i];
			++i;
		}
		
		if ( local_sum < (sum_of_elements - length + i) )// i == s - 1
		{
			comp[i]     -= 1;
			comp[i - 1] += 1;
		}
		else
		{
			--i;
			comp.back()  = comp[i] - 1;
			comp[i]      = 1;
			comp[i - 1] += 1;
		}
		
		return comp[0] == sum_of_elements - length + 1;
	}
	

	Compositions
	generate_compositions(unsigned sum_of_elements, unsigned length)
	{
		Compositions result;
		Composition  compos(length, 1);
		compos.back() = static_cast<int>(sum_of_elements - length + 1);
		
		while ( compos[0] != (sum_of_elements - length + 1) )
		{
			result.push_back(compos);
			get_next_composition(compos, sum_of_elements);
		}
		result.push_back(compos);
		return result;
		
	}
	
	std::vector<Dimensions>
	binomial_coefficient(size_t from, size_t to, size_t length)
	{
		size_t n = to - from + 1;
		std::vector<bool> sep(n, 0);
		std::fill(sep.begin(), sep.begin() + length, 1);
		
		std::vector<Dimensions> res;
		
		do
		{
			Dimensions dims;
			for (size_t i = 0; i < n; i++)
			{
				if (sep[i])
				{
					dims.push_back(i + from);
				}
			}
			res.push_back(dims);
		} while (std::prev_permutation(sep.begin(), sep.end()));
		
		return res;
	}
	
	bool
	unified_checker(int lip, int s, Callback const &composition_checker)
	{
		bool allright = true;
		int  dim_length = 2;
		
		while ( allright && dim_length <= s && dim_length <= lip )
		{
			auto dim_combs = binomial_coefficient(0, s - 1, dim_length);
			auto d_vecs    = generate_compositions(lip, dim_length);
			
			for (auto iter_dim_combs = dim_combs.begin(); allright && iter_dim_combs != dim_combs.end(); ++iter_dim_combs)
			{
				for (auto iter_d_vecs = d_vecs.begin(); allright && iter_d_vecs != d_vecs.end();)
				{
					allright = composition_checker(*iter_dim_combs, d_vecs, iter_d_vecs);
				}
			}
			
			++dim_length;
		}
		
		return allright;
	}
	
	tms::BasicInt
	impl_find_sobol_defect(int cur_t, int m, int s, Callback const &composition_checker)
	{
		while ( cur_t > 0 && unified_checker(m - cur_t, s, composition_checker) )
		{
			--cur_t;
		}
		
		return cur_t + 1;
	}
	
};



tms::BasicInt
tms::analysis::find_sobol_defect(tms::Sobol const &net, tms::analysis::SobolTestType test_type, int t_estimate)
{
	unsigned t_star = net.t_estimate();
	
	t_estimate = ( t_estimate < 0 || t_estimate > t_star ) ? t_star : t_estimate;
	
	Callback callback_wrapper;
	
	if ( test_type == SobolTestType::poly )
	{
		ClosurePOLY callback(net);
		callback_wrapper = [&](Dimensions const &dims, Compositions const &d_vecs, Compositions::iterator &iter_d_vecs) {
			return callback(dims, d_vecs, iter_d_vecs);
		};
		return impl_find_sobol_defect(t_estimate + 1, net.m(), net.s(), callback_wrapper);
	}
	else if ( test_type == SobolTestType::lip )
	{
		ClosureLIP callback(net, 0);
		callback_wrapper = [&](Dimensions const &dims, Compositions const &d_vecs, Compositions::iterator &iter_d_vecs) {
			return callback(dims, d_vecs, iter_d_vecs);
		};
		return impl_find_sobol_defect(t_estimate + 1, net.m(), net.s(), callback_wrapper);
	}
	else if ( test_type == SobolTestType::lip_mod )
	{
		ClosureLIP callback(net, 1);
		callback_wrapper = [&](Dimensions const &dims, Compositions const &d_vecs, Compositions::iterator &iter_d_vecs) {
			return callback(dims, d_vecs, iter_d_vecs);
		};
		return impl_find_sobol_defect(t_estimate + 1, net.m(), net.s(), callback_wrapper);
	}
		
	return ~(0);
}
