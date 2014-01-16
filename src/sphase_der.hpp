#ifndef __SPHASE_HPP__
#define __SPHASE_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>

void selective_phase( Sequence::gsl_uniform & uni,
		      Sequence::gsl_uniform01 & uni01,
		      std::vector<Sequence::chromosome> & sample,
		      Sequence::arg & sample_history,
		      const int & ttl_nsam,
		      int * NSAM,
		      int * nlinks,
		      double * t,
		      const double & rho,
		      const int & nsites,
		      const int & N,
		      const double & f,
		      const double & s,
		      const std::vector<double> & path,
		      const int & X);

void selective_phase_competing( Sequence::gsl_uniform & uni,
				Sequence::gsl_uniform01 & uni01,
				Sequence::gsl_exponential & expo,
				std::vector<Sequence::chromosome> & sample,
				Sequence::arg & sample_history,
				const int & ttl_nsam,
				int * NSAM,
				int * nlinks,
				double * t,
				const double & rho,
				const int & nsites,
				const int & N,
				const double & f,
				const double & s,
				const std::vector<double> & path,
				const int & X);

#endif
