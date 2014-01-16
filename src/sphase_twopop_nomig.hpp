#ifndef __SPHASE_TWOPOP_NOMIG_HPP__
#define __SPHASE_TWOPOP_NOMIG_HPP__

#include<Sequence/Coalescent/SimTypes.hpp>
#include<Sequence/RNG/gsl_rng_wrappers.hpp>

void selective_phase_twopop_nomig( Sequence::gsl_uniform & uni,
				   Sequence::gsl_uniform01 & uni01,
				   std::vector<Sequence::chromosome> & sample,
				   Sequence::arg & sample_history,
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
