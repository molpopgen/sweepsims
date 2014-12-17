#ifndef __SPHASE_TWOPOP_NOMIG_HPP__
#define __SPHASE_TWOPOP_NOMIG_HPP__

#include <common.hpp>
#include<Sequence/Coalescent/SimTypes.hpp>
#include<functional>

void selective_phase_twopop_nomig( std::function<double(const double &,const double &)> & uni,
				   std::function<double()> & uni01,
				   std::vector<chromosome> & sample,
				   ARG & sample_history,
				   const int & ttl_nsam,
				   int * NSAM,
				   int * nlinks,
				   double * t,
				   bool * merged,
				   const double & rho,
				   const int & nsites,
				   const int & N,
				   const double & tmerge,
				   const double & s,
				   const std::vector<double> & path,
				   const int & X,
				   const int & k);
#endif
