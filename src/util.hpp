#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <common.hpp>
#include <Sequence/Coalescent/SimTypes.hpp>
#include <cassert>

int slinks( const chromosome & c,
	    const int & X );

std::pair<int,int> 
pick_uniform_spot2(const double & random_01,
		   const int & nlinks,
		   const int & X,
		   std::vector<chromosome>::const_iterator sample_begin,
		   const unsigned & current_nsam,
		   const int & deme);

std::pair<int,int> 
pick_uniform_spot3(const double & random_01,
		   const int & nlinks,
		   std::vector<chromosome>::const_iterator sample_begin,
		   const unsigned & current_nsam,
		   const int & deme);
#endif
