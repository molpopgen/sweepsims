#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/PolySIM.hpp>
#include <detpath.hpp>
#include <sphase.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

int main( int argc, char **argv )
{
  if(argc != 12)
    {
      cerr << "sweep_stochCG_soft : single soft sweep at time tau. Stochastic trajectories using method of Coop and Griffiths\n"
	   << "usage: sweep_stochCG nsam nreps N 2Ns tau f0 X rho_site nsites theta_site seed\n"
	   << "where:\n"
	   << "\ttau is in units of 4N generations\n"
	   << "\tX is the position of the selected site, from 1 to nsites, inclusive.  if X is -1, it is random from 1 to nsites for each replicate\n";
      exit(1);
    }

  const int nsam = atoi(argv[1]);
  const int nreps = atoi(argv[2]);
  const int N = atoi(argv[3]);
  const double alpha = atof(argv[4]);
  const double tau = atof(argv[5]);
  const double f0=atof(argv[6]);
  const int X = atoi(argv[7])-1; //user should input 1 to nsites
  const double rho = atof(argv[8]);
  const int nsites = atoi(argv[9]);
  const double theta = atof(argv[10]);
  const double s = alpha/(2.*N);
  const unsigned seed = atoi(argv[11]);

  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';

  vector<chromosome> initialized_sample = init_sample( vector<int>(1,nsam), 
						       nsites );
  marginal initialized_marginal = init_marginal(nsam);
  vector<double> path,path2;
  path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  double ar;

  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };

  double rcoal,rrec,tcoal,trec,t;
  const int k = 10;
  const double dtp = 1./double(2*k*N),ifreq=f0;
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      int NSAM = nsam;
      int nlinks = NSAM*(nsites-1);
      vector<chromosome> sample(initialized_sample);
      ARG sample_history(1,initialized_marginal);
      bool neutral = true;
      t=0.;
      while ( NSAM > 1 )
	{
	  if( neutral )
	    {
	      rcoal = double(NSAM*(NSAM-1));
	      rrec = rho*double(nlinks);
	      
	      tcoal = expo(1./rcoal);
	      trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	      double tmin = min(tcoal,trec);
	      if ( (t==0.&&tau==0.) || (t < tau && t+tmin >= tau) )
		{
		  t = tau;
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
	      int sel_site = (X>=0) ? X : int(gsl_ran_flat(r,0.,nsites));
	      cerr  << sel_site << '\n';
	      ConditionalTraj(uni01,&path,N,s,dtp,ifreq);
	      ConditionalTrajNeutral(uni01,&path2,dtp,f0);
	      /*
	      for(unsigned i=0;i<path.size();++i)
		{
		  cout << path[i] << '\n';
		}
	      for(unsigned i=0;i<path2.size();++i)
		{
		  cout << path2[i] << '\n';
		}
	      exit(1);
	      */
	      copy(path2.begin(),path2.end(),
		   back_inserter(path));
	      const double told = t;
	      selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,
			      &nlinks,&t,rho,nsites,N,path,sel_site,1./double(2*N),1./double(4*k*N),k);//1./double(2*N));
	      neutral = true;
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,sample_history,theta*double(nsites));
      cout << d << '\n' ;
    }
}
