#ifndef __SWEEP_DER_ARG_HPP__
#define __SWEEP_DER_ARG_HPP__
#include <Sequence/Coalescent/SimTypes.hpp>
#include <gsl/gsl_rng.h>

/*!
  Simulate a selective sweep in a derived population that has undergone
  a bottleneck.
  \param r a gsl rng
  \param n1 sample size of  pop
  \param nsites number of sites
  \param rho 4Nr/site
  \param tr time of recovery from bottleneck, units of 4N0 generations (units of current ancestral pop size).
  \param d duration of bneck, same units as tr
  \param f Nb/N0
  \param N take this to mean N0
  \param s selection coefficient, such that alpha = 2*N0*f*s in the derived pop
  \param tau time of sweep, same units as tr
  \param X position of beneficial mutation
  \param path trajectory of beneficial allele, calculated with xi=1/(2N0f) and alpha = 2N0fs
*/
Sequence::arg sweep_bneck_arg(gsl_rng * r,
			      const int & n1,
			      const int & nsites,
			      const double & rho,
			      const double & tr,
			      const double & d,
			      const double & f,
			      const int & N,
			      const double & s,
			      const double & tau,
			      const int & X,
			      const std::vector<double> & path,
			      const std::vector<Sequence::chromosome> & initialized_sample,
			      const Sequence::marginal & initialized_marginal,
			      const double & dt);

Sequence::arg sweep_bneck_argCG(gsl_rng * r,
				const int & n1,
				const int & nsites,
				const double & rho,
				const double & tr,
				const double & d,
				const double & f,
				const int & N,
				const double & s,
				const double & tau,
				const int & X,
				const std::vector<double> & path,
				const std::vector<Sequence::chromosome> & initialized_sample,
				const Sequence::marginal & initialized_marginal,
				const int & k,
				const double & dt);

#endif
