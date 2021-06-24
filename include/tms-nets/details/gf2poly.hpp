/**
 * @file    gf2poly.hpp
 *
 * @author Vadim Piven
 * @author Alexey Burimov
 */

#ifndef TMS_NETS_GF2POLY_HPP
#define TMS_NETS_GF2POLY_HPP

#include "common.hpp"

#include <map>

/// @todo exceptions

// coefficient number of polynomial over GF(2) is a such integer number that it's n-th bit is equal to n-th coefficient of polynomial.

/** @namespace tms::gf2poly
 *  @brief Contains specific polynomial-related functions that are necessary for digital nets construction */
namespace tms::gf2poly
{
	/** Returns polynomial over GF(2) with the specified coefficients (represents a wrapper above irrpoly::gfpoly constructor).
	 *  @param [in] coeffs - desired polynomial coefficients */
	Polynomial              make_gf2poly(std::vector<uintmax_t> const &coeffs);
	/** Returns a polynomial which i-th coefficient is equal to i-th bit of a given number
	 *  @param [in] coeffs_number - number corresponding to  */
	Polynomial              number_to_poly(uintmax_t coeffs_number);
	/** Returns a number which i-th bit is equal to i-th coefficient of a given polynomial
	 *  @param [in] coeffs_number - number corresponding to  */
	uintmax_t               poly_to_number(Polynomial const &poly);
	/** Returns the value of polynomial raised to the power
	 *  @param [in] poly – base polynomial
	 *  @param [in] expo – integer power value */
	Polynomial              pow(irrpoly::gfpoly poly, uintmax_t expo);
	/** Returns the of polynomial raised to the power
	 *  */
	Polynomial              find_inverse_straightforward(Polynomial arg, Polynomial const &mod);
	/** Returns the of polynomial raised to the power
	 *  */
	Polynomial              find_inverse(Polynomial arg, Polynomial mod);
	/** Generates vector of first least-degree irreducible polynomials over GF(2).
	 *  @param [in] amount - amount of irreducible polynomials to generate
	 *  @param [in] degree_limit - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys(unsigned int amount, unsigned int degree_limit = max_nbits);
	/** Generates vector of first least-degree irreducible polynomials over GF(2) using multithreading.
	 *  @param [in] amount - amount of irreducible polynomials to generate
	 *  @param [in] degree_limit - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys_in_parallel(unsigned int amount, unsigned int degree_limit = max_nbits);
	/** Generates vector of first irreducible polynomials over GF(2) with specified degrees.
	 *  @param [in] degrees - vector of desired degrees of irreducible polynomials
	 *  @param [in] degree_limit - upper limit for sum of degrees of polynomials */
	std::vector<Polynomial> generate_irrpolys_with_degrees(std::vector<unsigned int> const &degrees, unsigned int degree_limit = max_nbits);
	
	/** Generates vector of all first irreducible polynomials over GF(2) with degrees <= degree in "ascending" order
	 *  @param [in] degree - greatest degree of generated polynomials */
	std::vector<Polynomial> generate_all_irrpolys_until_degree(unsigned int degree);
	
	std::vector<Polynomial> generate_all_irrpolys_of_degree(unsigned int degree);
	
};


#endif /* gf2poly_hpp */
