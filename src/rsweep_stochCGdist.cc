#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/RNG/gsl_rng_wrappers.hpp>
#include <Sequence/PolySIM.hpp>
#include <boost/bind.hpp>
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
  if(argc != 9)
    {
      cerr << "rsweep_stochCG: recurrent hitch-hiking within the sampled region, and at linked sites (following Jensen et al. 2007, Genetics), with trajectories determined by the method of Coop and Griffiths\n"
	   << "usage: rsweep_linked_stoch nsam nreps N mean_s lambda rho_locus nsites theta_locus seed\n"
	   << "where:\n"
	   << "\tmean_s is used to choose s for each sweep.  This is done by sampling a value of s from an exponential with mean mean_s, and then accepting it with probability 2s.  A large number of acceptances would be extreme-value distributed\n"
	   << "\tlambda = expected # of sweeps per site per 4N generations\n";
      exit(1);
    }

  int argn=1;
  const int nsam = atoi(argv[argn++]);
  const int nreps = atoi(argv[argn++]);
  const int N = atoi(argv[argn++]);
  const double mean_s = atof(argv[argn++]);
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


  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  gsl_uniform01 uni01(r); 
  gsl_uniform uni(r);     
  gsl_exponential expo(r);
  gsl_poisson poiss(r); 

  double rcoal,rrec,rsweepin=SEQMAXDOUBLE,rsweepout=SEQMAXDOUBLE,tcoal,trec,t,tsweep;
  const int k = 10;
  const double dtp = 1./double(2*k*N),ifreq=1./double(2*N);
  double rdm;
  std::cerr << "//\n";
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      int NSAM = nsam;
      int nlinks = NSAM*(nsites-1);
      vector<chromosome> sample(initialized_sample);
      Sequence::arg sample_history(1,initialized_marginal);
      bool neutral = true;
      t=0.;
      double s = gsl_ran_exponential(r,mean_s);
      rdm = gsl_rng_uniform(r);
      while( rdm > 2.*s )
	{
	  s = gsl_ran_exponential(r,mean_s);
	  rdm = gsl_rng_uniform(r);
	}
      path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));
      while ( NSAM > 1 )
	{
	  double maxd = 4.*double(N)*s;
	  if( neutral )
	    {
	      rcoal = double(NSAM*(NSAM-1));
	      rrec = rho*double(nlinks);

	      //rate of sweeps within locus
	      rsweepin = lambda*double(nsites);
	      //rate of sweeps at linked sites
	      rsweepout = 2.*lambda*maxd/rho;
	      tcoal = expo(1./rcoal);
	      trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	      tsweep = (rsweepin+rsweepout>0.)?expo(1./(rsweepin+rsweepout)): SEQMAXDOUBLE;
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
	      std::cerr << t << '\t' << NSAM << '\t';
	      int X = 0;
	      //is sweep within region, or at nearby linked sites?
	      if( gsl_rng_uniform(r) < rsweepin/(rsweepin+rsweepout) )
		{
		  X = int(gsl_ran_flat(r,0.,nsites));
		}
	      else
		{
		  //X is now distance from neutral region to selected site
		  X = int((gsl_rng_uniform(r)*maxd)/rho);
		  
		  //is X to left or right of region?
		  if( gsl_rng_uniform(r) < 0.5 )
		    {
		      X *= -1; // to the left
		    }
		  else
		    {
		      X += nsites; // to the right
		    }
		}
	      std::cerr << X << '\t' << s << endl;
	      ConditionalTraj(boost::bind(gsl_rng_uniform,r),&path,N,s,dtp,ifreq);
	      selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,
			      &nlinks,&t,rho,nsites,N,path,X,1./double(2*N),1./double(4*k*N));
	      neutral = true;
	      //set s for the next sweep
	      s = gsl_ran_exponential(r,mean_s);
	      while( rdm > 2.*s )
		{
		  s = gsl_ran_exponential(r,mean_s);
		  rdm = gsl_rng_uniform(r);
		}

	      path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,
					  sample_history,theta);
      cout << d << '\n' ;
    }
}





