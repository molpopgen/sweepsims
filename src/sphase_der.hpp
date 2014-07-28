#ifndef __SPHASE_HPP__
#define __SPHASE_HPP__
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <functional>
#include <Sequence/Coalescent/SimTypes.hpp>

void selective_phase( std::function<double(const double&,const double&)> & uni,
		      std::function<double()> & uni01, 
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

void selective_phase_competing( std::function<double(const double&,const double&)> & uni,
				std::function<double()> & uni01,
				std::function<double(const double &)> & expo,
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
