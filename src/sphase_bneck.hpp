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
		      const double & f,
		      const double & rho,
		      const int & nsites,
		      const int & N,
		      const double & s,
		      const std::vector<double> & path,
		      const int & X,
		      const double & dt);

void selective_phaseCG( Sequence::gsl_uniform & uni,
			Sequence::gsl_uniform01 & uni01,
			std::vector<Sequence::chromosome> & sample,
			Sequence::arg & sample_history,
			const int & ttl_nsam,
			int * NSAM,
			int * nlinks,
			double * t,
			const double & f,
			const double & rho,
			const int & nsites,
			const int & N,
			const double & s,
			const std::vector<double> & path,
			const int & X,
			const int & k,
			const double & dt);

#endif
