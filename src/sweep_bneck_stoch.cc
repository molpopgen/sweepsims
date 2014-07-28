#include <Sequence/Coalescent/Initialize.hpp>
#include <Sequence/Coalescent/TreeOperations.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/PolySIM.hpp>
#include <detpath.hpp>
#include <sweep_bneck_arg.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>
#include <cmath>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

int main( int argc, char **argv )
{
  if (argc != 14 )
    {
      cerr << "sweep_bneck: sweep during a bottleneck. Stochastic trajectories using forward simulation\n"
	   << "usage: sweep_bneck n1 nsam nreps t_r d f N 2Ns tau X rho_site nsites theta_site seed\n"
	   << "where\n"
	   << "\tt_r = time of recovery from bneck, units of 4N0 (N0 = ancestral effective size)\n"
	   << "\td = duration of bneck, units of 4N0\n"
	   << "\tf = severity of bneck.  Size of derived pop. during bneck is fN0\n"
	   << "\ttau = time of sweep, units of 4N0\n"
	   << "\tX is position of selected site, from 1 to nsites\n";
      exit(1);
    }
  const int n1 = atoi(argv[1]);
  const int nsam = n1;
  const int nreps = atoi(argv[2]);
  const double tr = atof(argv[3]);
  const double d = atof(argv[4]);
  const double f = atof(argv[5]);
  const int N = atoi(argv[6]);
  const double alpha = atof(argv[7]);
  const double tau = atof(argv[8]);
  const int X = atoi(argv[9])-1; //user should input 1 to nsites
  const double rho = atof(argv[10]);
  const int nsites = atoi(argv[11]);
  const double theta = atof(argv[12]);
  const unsigned seed = atoi(argv[13]);
  const double s = alpha/(2.*double(N)*f);
  //check that sweep params are sane
  if( tau < tr || tau >= tr+d )
    {
      cerr << "time of sweep not valid: ";
      if(tau<tr)
	{
	  cerr << "time of sweep < recovery time from bneck (tr)\n";
	}
      else
	{
	  cerr << "time of sweep is after populations have merged\n";
	}
      exit(1);
    }
  if( tau + (f*-std::log(1./double(2*f*N))/double(2*int(f*N)*s)) >= (tr+d) )
    {
      cerr << "sweep occurs during bottleneck but will "
	   << "end after bottleneck begins\n";
      exit(1);
    }
  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';
  vector<int> c(1,n1);
  vector<chromosome> initialized_sample = init_sample( c,nsites );
  marginal initialized_marginal = init_marginal(nsam);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };

  double ar;
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      int sel_site = X;
      if(X < 0)
	{
	  sel_site = int(gsl_ran_flat(r,0.,double(nsites)));
	}
      vector<double> path;
      do
	{
	  path = rev_traj(r,1.-1./(2.*f*double(N)),-2.*s,0.5,int(f*double(N)),&ar);
	}
      while( (tr+(f*double(path.size())/(4.*f*double(N)))) >= tr+d );
      Sequence::arg sample_history =  sweep_bneck_arg(r, n1, nsites, rho,
						      tr, d, f, N, s, tau, sel_site, path,
						      initialized_sample,initialized_marginal,
						      1./floor(4.*f*double(N)));
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,sample_history,theta*double(nsites));
      cout << d << endl ;
      if(X<0) cerr << sel_site << '\n';
    }
}






