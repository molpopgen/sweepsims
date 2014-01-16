#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/RNG/gsl_rng_wrappers.hpp>
#include <Sequence/PolySIM.hpp>
#include <boost/bind.hpp>
//#include <detpath.hpp>
#include <sphase.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>

using namespace std;
using namespace Sequence;

int main( int argc, char **argv )
{
  if( argc != 10 )
    {
      cerr << "rsweep_within_stochCG: recurrent HH, only within the sampled region (!!!!!).  Stochastic trajectories using the method of Coopy and Griffiths\n"
	   << "usage: rsweep_within_stochCG nsam nreps N s lambda rho_locus nsites theta seed\n"
	   << "where:\n"
	   << "lambda is the 1/(mean time to next sweep).  Note: very different from Braverman et al., Przeworski, and the rest of the world.  You have been warned\n";
      exit(1);
    }
  int argn=1;
  const int nsam = atoi(argv[argn++]);
  const int nreps = atoi(argv[argn++]);
  const int N = atoi(argv[argn++]);
  const double s = atof(argv[argn++]);
  const double lambda = atof(argv[argn++]);
  const double rho_locus = atof(argv[argn++]);
  const int nsites = atoi(argv[argn++]);
  const double theta = atof(argv[argn++]);
  const unsigned seed = atoi(argv[argn++]);
  const double rho = rho_locus/double(nsites-1);
  if(lambda <= 0.) 
    {
      cerr << "error: lambda must be > 0!\n";
      exit(1);
    }
  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';
  vector<chromosome> initialized_sample = init_sample( vector<int>(1,nsam), 
						       nsites );
  marginal initialized_marginal = init_marginal(nsam);
  vector<double> path;
  path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  gsl_uniform01 uni01(r); 
  gsl_uniform uni(r);     
  gsl_exponential expo(r);
  gsl_poisson poiss(r); 

  double rcoal,rrec,rsweepin=SEQMAXDOUBLE,tcoal,trec,t,tsweep;
  const int k = 2;
  const double dtp = 1./double(2*k*N),ifreq=1./double(2*N);
  //unsigned nsweeps = 0;
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      int NSAM = nsam;
      int nlinks = NSAM*(nsites-1);
      vector<chromosome> sample(initialized_sample);
      Sequence::arg sample_history(1,initialized_marginal);
      bool neutral = true;
      t=0.;
      while ( NSAM > 1 )
	{
	  if( neutral )
	    {
	      rcoal = double(NSAM*(NSAM-1));
	      rrec = rho*double(nlinks);

	      //rate of sweeps within locus
	      rsweepin = lambda*double(nsites);
	      tcoal = expo(1./rcoal);
	      trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	      tsweep = (rsweepin>0.)?expo(1./(rsweepin)): SEQMAXDOUBLE;
	      double tmin = min(tcoal,trec);
	      if( tsweep <= tmin )
		{
		  t += tsweep;
		  neutral = false;
		}
	      else // neutral phase events
		{
		  t += tmin;
		  if ( tcoal < trec ) //coalescence
		    {
		      //std::pair<int,int> two = pick2(uni,NSAM);
		      std::pair<int,int> two = pick2(uni,NSAM);
		      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
				       &nlinks,&sample,&sample_history);
		    }
		  else //recombination
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(gsl_rng_uniform(r),
								     nlinks,
								     sample.begin(),NSAM);
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		    }
		  if (NSAM < int(sample.size())/5)
		    {
		      sample.erase(sample.begin()+NSAM+1,sample.end());
		    }
		}
	    }
	  else //selective phase
	    {
	      int X = int(gsl_ran_flat(r,0.,nsites));
	      ConditionalTraj(uni01,&path,N,s,dtp,ifreq);
	      selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,
			      &nlinks,&t,rho,nsites,N,path,X,1./double(2*N),1./double(4*k*N),k);
	      neutral = true;
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,
					  sample_history,theta);
      cout << d << '\n' ;
    }
}





