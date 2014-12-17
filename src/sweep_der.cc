#include <Sequence/Coalescent/Initialize.hpp>
#include <Sequence/Coalescent/TreeOperations.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/PolySIM.hpp>
#include <detpath.hpp>
#include <sweep_der_arg.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>
#include <cmath>
#include <fstream>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;
using namespace Sequence;

// struct params
// {
//   int n1,n2,nreps,N,X,nsites;
//   double tr,d,f,alpha,tau,rho,theta;
//   params(const int & _n1,const int & _n2,
// 	 const int & _nreps, const int & _N,
// 	 const int & _X, const int _nsites,
// 	 const double & _tr, const double & _d,
// 	 const double & _f, const double & _alpha, const double _tau,
// 	 const double & _rho, const double _theta)
//     : n1(_n1),n2(_n2),nreps(_nreps),N(_N),X(_X),nsites(_nsites),
//       tr(_tr),d(_d),f(_f),alpha(_alpha),tau(_tau),
//       rho(_rho),theta(_theta)
//   {
//   }
// };

//params read_params( ifstream & in );

struct sortpos : public binary_function<polymorphicSite,polymorphicSite,bool>
{
  inline bool operator()(const polymorphicSite & a, const polymorphicSite & b) const
  {
    return a.first < b.first;
  }
};

int main( int argc, char **argv )
{
  if (argc != 15)
    {
      cerr << "sweep_der: sweeps in a derived population, following the 2-population model of Thornton and Jensen (2007) Genetics.  Trajectories are deterministic\n"
	   << "usage: sweep_der nsam_ancestral nsam_derived nreps t_r d f N 2Ns tau X rho_site nsites theta_site seed\n"
	   << "where\n"
	   << "\tt_r = time of recovery from bneck, units of 4N0 (N0 = ancestral effective size)\n"
	   << "\td = duration of bneck, units of 4N0\n"
	   << "\tf = severity of bneck.  Size of derived pop. during bneck is fN0\n"
	   << "\ttau = time of sweep, units of 4N0\n"
	   << "\tX is the position of the selected site, from 1 to nsites, inclusive\n"
	   << "note: the fixed difference that results from the selected site is not in the output!\n";
      exit(1);
    }
  const int n1 = atoi(argv[1]);
  const int n2 = atoi(argv[2]);
  const int nsam = n1+n2;
  const int nreps = atoi(argv[3]);
  const double tr = atof(argv[4]);
  const double d = atof(argv[5]);
  const double f = atof(argv[6]);
  const int N = atoi(argv[7]);
  const double alpha = atof(argv[8]);
  const double tau = atof(argv[9]);
  const int X = atoi(argv[10])-1; //user should input 1 to nsites
  const double rho = atof(argv[11]);
  const int nsites = atoi(argv[12]);
  const double theta = atof(argv[13]);
  const unsigned seed = atoi(argv[14]);
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
	   << "end after populations merge\n";
      exit(1);
    }

  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';
  vector<int> c(1,n1);
  c.push_back(n2);
  vector<chromosome> initialized_sample = init_sample( c,nsites );
  marginal initialized_marginal = init_marginal(nsam);
  vector<double> path = detpath(int(N*f),s);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  //gsl_rng_set(r,time(0));
  gsl_rng_set(r,seed);

  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };

  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      cerr << rep << '\n';
      ARG sample_history =  sweep_der_arg(r, n1, n2, nsites, rho,
						    tr, d, f, N, s, tau, X, path,
						    initialized_sample,initialized_marginal);
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,sample_history,theta*double(nsites));

      //The infinite_sites function doesn't add in a fixed difference corresponding to the selected site
      //so, we do that here:
      vector<polymorphicSite> vps(d.sbegin(),d.send());
      vps.push_back( polymorphicSite( double(X)/double(nsites), string(n1,'0')+string(n2,'1') ) );
      sort(vps.begin(),vps.end(),sortpos());
      
      d.assign(vps.begin(),vps.end());
      cout << d << endl ;
    }
}






