#ifndef __SPHASE_HPP__
#define __SPHASE_HPP__
#include <common.hpp>
#include <Sequence/Coalescent/SimTypes.hpp>
#include <functional>
void selective_phase( std::function<double(const double &,const double &)> & uni,
		      std::function<double()> & uni01,
		     std::vector<chromosome> & sample,
		     ARG & sample_history,
		     const int & ttl_nsam,
		     int * NSAM,
		     int * nlinks,
		     double * t,
		     const double & rho,
		     const int & nsites,
		     const int & N,
		     const std::vector<double> & path,
		     const int & X,
		     const double & pmin,
		     const double & dt,
		     const int & k = 1);

void selective_phase_partial(std::function<double(const double &,const double &)> & uni,
			     std::function<double()> & uni01,
			     std::vector<chromosome> & sample,
			     ARG & sample_history,
			     const int & ttl_nsam,
			     int * NSAM,
			     int * nlinks,
			     double * t,
			     const double & rho,
			     const int & nsites,
			     const int & N,
			     const std::vector<double> & path,
			     const int & X,
			     const double & starting_freq,
			     const double & pmin,
			     const double & dt,
			     const int & k = 1);

void selective_phase_competing(std::function<double(const double &,const double &)> & uni,
			       std::function<double()> & uni01,
			       std::function<double(const double &)> & expo,
			       std::vector<chromosome> & sample,
			       ARG & sample_history,
			       const int & ttl_nsam,
			       int * NSAM,
			       int * nlinks,
			       double * t,
			       const double & rho,
			       const int & nsites,
			       const int & N,
			       const std::vector<double> & path,
			       const int & X,
			       const double & pmin,
			       const double & dt,
			       const int & k = 1);

/*
void selective_phase_linked(gsl_uniform & uni,
			    gsl_uniform01 & uni01,
			    std::vector<chromosome> & sample,
			    ARG & sample_history,
			    const int & ttl_nsam,
			    int * NSAM,
			    int * nlinks,
			    double * t,
			    const double & rho_sel,
			    const double & rho,
			    const int & nsites,
			    const std::vector<double> & path,
			    const double & dt,
			    const double & pmin,
			    const bool & isleft,
			    const int & k = 1);
*/
/*
void selective_phase_CG( gsl_uniform & uni,
			 gsl_uniform01 & uni01,
			 std::vector<chromosome> & sample,
			 ARG & sample_history,
			 const int & ttl_nsam,
			 int * NSAM,
			 int * nlinks,
			 double * t,
			 const double & rho,
			 const int & nsites,
			 const int & N,
			 const std::vector<double> & path,
			 const double & pathdt,
			 const int & k,
			 const int & X);
*/
#endif
