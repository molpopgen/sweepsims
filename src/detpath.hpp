#ifndef __DETPATH_HPP__
#define __DETPATH_HPP__

#include <vector>
#include <gsl/gsl_rng.h>

std::vector<double> detpath(const unsigned & N,
			    const double & s,
			    const unsigned N_for_delta_p = 0);

std::vector<double> rev_traj( gsl_rng * r,
			      const double & initial_freq,
			      const double & s,
			      const double & h,
			      const unsigned & N,
			      double * acc_rate);

std::vector<double> rev_traj_haploid( gsl_rng * r,
				     const double & initial_freq,
				     const double & s,
				     const unsigned & N);
#endif
