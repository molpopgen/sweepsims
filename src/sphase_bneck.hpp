#ifndef __SPHASE_HPP__
#define __SPHASE_HPP__

#include <Sequence/Coalescent/SimTypes.hpp>
#include <functional>
void selective_phase (std::function<double(const double &, const double &)> & uni,
		      std::function<double()> & uni01,
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

void selective_phaseCG(std::function<double(const double &, const double &)> & uni,
		       std::function<double()> & uni01,
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
