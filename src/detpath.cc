#include <detpath.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <cstdlib>
std::vector<double> detpath(const unsigned & N,
			    const double & s,
			    const unsigned N_for_delta_p)
{
  double p= (N_for_delta_p == 0) ? 1./double(2*N) : 1./double(2*N_for_delta_p);
  int m = int(std::floor(4*std::log(2.*double(N))/s));
  if(m == 0)
    {
      std::cerr << "fatal error: trajectory is of length 0\n";
      exit(10);
    }
  std::vector<double> path(m,0.);
  path[0]=p; 
  int ngen=1;
  while (path[ngen] < 1) {
    path[ngen] = (double) p/(p + (1-p)*std::exp(-s*ngen));
    if (path[ngen] > (1-p))  path[ngen]=1.0;
    else ngen++; 
  }
  //std::reverse(path.begin(),path.begin()+ngen+1);
  //path.erase(path.begin()+ngen+1,path.end());
  std::reverse(path.begin(),path.begin()+ngen+1);
  path.erase(path.begin()+ngen+1,path.end());
  return path;
}

std::vector<double> rev_traj( gsl_rng * r,
			      const double & initial_freq,
			      const double & s,
			      const double & h,
			      const unsigned & N,
			      double * acc_rate)
{
  const double w11 = 1. + s, w10 = 1. + s*h, w00=1.;
  double p1,mf;
  bool accepted = false;
  unsigned ntries=0;
  while (1)
    {
      ++ntries;
      std::vector<double> traj(1,initial_freq);
      double x = initial_freq;
      while( std::fabs( 1. - x ) > std::numeric_limits<double>::epsilon()
	     && std::fabs( 0. - x ) > std::numeric_limits<double>::epsilon() )
	{
	  mf = (w11*x*x) + (w10*2.*x*(1.-x)) + (w00*(1.-x)*(1.-x));
	  p1 = x*(w11*x + w10*(1.-x))/mf;
	  //std::cerr << x << ' ' << p1 << '\n';
	  unsigned n = gsl_ran_binomial(r,p1,2*N);
	  x  = double(n)/double(2*N);
	  traj.push_back(x);
	}
      if(  std::fabs( 0. - x ) < std::numeric_limits<double>::epsilon() )
	{
	  *acc_rate = 1./double(ntries);
	  return traj;
	  //accepted = true;
	}
    }
  *acc_rate = std::strtod("NAN",NULL);
  return std::vector<double>();
}

std::vector<double> rev_traj_haploid( gsl_rng * r,
				      const double & initial_freq,
				      const double & s,
				      const unsigned & N)
{
  const double w1 = 1. + s, w0 = 1.;
  double p1,mf;
  bool accepted = false;
  unsigned ntries=0;
  double x = initial_freq;
  std::vector<double> traj;
  do
    {
      ++ntries;
      traj = std::vector<double>(1,initial_freq);
      x = initial_freq;
      while( std::fabs( 1. - x ) > std::numeric_limits<double>::epsilon()
	     && std::fabs( 0. - x ) > std::numeric_limits<double>::epsilon() )
	{
	  mf = (w1*x) + (w0*(1.-x));
	  p1 = (w1*x)/mf;
	  unsigned n = gsl_ran_binomial(r,p1,N);
	  x  = double(n)/double(N);
	  traj.push_back(x);
	}
    }
  while( std::fabs( 1. - x ) > std::numeric_limits<double>::epsilon() );
  //  *acc_rate = 1./double(ntries);
  std::reverse(traj.begin(),traj.end());
  return traj;
}
