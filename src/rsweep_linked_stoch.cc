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
  if(argc != 10)
    {
      cerr << "rsweep_linked_stoch: recurrent hitch-hiking at linked sites(i.e. Braverman 1995), with trajectories determined by forward simulation\n"
	   << "usage: rsweep_linked_stoch nsam nreps N s lambda rho_locus nsites theta_locus seed\n"
	   << "where:\n"
	   << "\tlambda = expected # of sweeps per site per 4N generations\n";
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
  const double rho = rho_locus/double(nsites-1);
  const double theta_locus = atof(argv[argn++]);
  //  const double s = alpha/(2.*double(N));
  const unsigned seed = atoi(argv[argn++]);

  std::copy(argv,argv+argc,std::ostream_iterator<char *>(cout," "));
  cout << endl;
  if(lambda <= 0.) 
    {
      cerr << "error: lambda must be > 0!\n";
      exit(1);
    }
  //const double rho = irho*double(nsites)/double(nsites-1);
  vector<chromosome> initialized_sample = init_sample( vector<int>(1,nsam), 
						       nsites );
  marginal initialized_marginal = init_marginal(nsam);
  double ar;
  //vector<double> path = detpath(N,s);
  /*
  copy(path.begin(),path.end(),
       ostream_iterator<double>(cout,"\n"));
  exit(1);
  */
  /*
  for(unsigned gen=0;gen<path.size();++gen)
    {
      cout << gen << ' ' <<path[gen]<<'\n';
    }
  exit(10);
  */
  //vector<double> genetic_map(nsites,rho);
  //genetic_map[nsites-1]=0.; //no rec b/w last link and stuff to the right

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); }; 

  double rcoal,rrec,rsel,tcoal,trec,t,tsweep,runi;
  double maxd = 4.*double(N)*s;
  bool left;
  int ssite;
  //double rho_locus = rho*double(nsites-1);
  //unsigned nsweeps = 0;
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
	      rsel = 2.*lambda*maxd/rho;
	      
	      trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	      tcoal = expo(1./rcoal);
	      tsweep = expo(1./rsel);
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
		      if (NSAM < int(sample.size())/5)
			{
			  sample.erase(sample.begin()+NSAM+1,sample.end());
			}
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
		} 
	    }
	  else //selective phase
	    {
	      //the fitness scheme is genic selection with fitnesses 1, 1+s, and 1+2s,
	      //which is equivalent to a model w/dominance and fitnesses 1, 1+hk, and 1+k, where k=2s and h=0.5.
	      left = (gsl_rng_uniform(r) < 0.5) ? true:false;
	      ssite = int((gsl_rng_uniform(r)*maxd)/rho);
	      if(left)
		{
		  ssite *= -1;
		}
	      else
		{
		  ssite += nsites;
		}
	      //cerr << ssite << '\n';
	      vector<double> path = rev_traj(r,1.-1./double(2*N),-2.*s,0.5,N,&ar);
	      /*
	      selective_phase_linked(uni,uni01,sample,sample_history,nsam,&NSAM,
				     &nlinks,&t,(gsl_rng_uniform(r)*maxd),
				     rho,nsites,path,1./double(4*N),1./double(2*N),left,1);
	      */
	      
	      
	      selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,&nlinks,&t,rho,nsites,N,path,
			      ssite,1./double(2*N),1./double(4.*N));
		
	      neutral = true;
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,sample_history,theta_locus);
      cout << d << endl; ;
    }
}





