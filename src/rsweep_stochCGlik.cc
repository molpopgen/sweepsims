#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/PolySIM.hpp>
#include <sphase.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

int main( int argc, char **argv )
{
  if(argc != 14)
    {
      cerr << "rsweep_stochCGcomplik: recurrent hitch-hiking within the sampled region, and at linked sites (following Jensen et al. 2007, Genetics), with trajectories determined by the method of Coop and Griffiths\n"
	   << "The selective phase is done using exponential distributions, rather than the Braverman algorithm\n"
	   << "usage: rsweep_linked_stoch nsam nreps N s lambda rho_locus nsites theta_locus obs_pi_site abs_tolerance rel_tolerance outfilename seed\n"
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
  const double theta = atof(argv[argn++]);
  const double obspi = atof(argv[argn++]);
  const double abs_tolerance = atof(argv[argn++]);
  const double rel_tolerance = atof(argv[argn++]);
  const char * ofn = argv[argn++];
  const unsigned seed = atoi(argv[argn++]);
  const double rho = rho_locus/double(nsites-1);
  if(lambda <= 0.) 
    {
      cerr << "error: lambda must be > 0!\n";
      exit(1);
    }
  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';

  ofstream ofs(ofn);
  if (! ofs)
    {
      cerr << "couldn't open " << ofn << '\n';
      exit(1);
    }
  ofs.close();
  vector<chromosome> initialized_sample = init_sample( vector<int>(1,nsam), 
						       nsites );
  marginal initialized_marginal = init_marginal(nsam);
  vector<double> path;
  path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  /*
  FILE * fp = fopen("seedstate","r");
  cerr << (fp == NULL) << '\n';
  if (!(fp == NULL))
    {
      cerr << "setting seed state\n";
      gsl_rng_fread(fp,r);
      fclose(fp);
    }
  */
  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); }; 

  double rcoal,rrec,rsweepin=SEQMAXDOUBLE,rsweepout=SEQMAXDOUBLE,tcoal,trec,t,tsweep;
  double maxd = 4.*double(N)*s;
  const int k = 10;
  const double dtp = 1./double(2*k*N),ifreq=1./double(2*N);
  unsigned prob = 0,prob2=0;
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      /*
      cerr << rep << endl;
      fp = fopen("seedstate","w");
      if (fp == NULL) abort();
      gsl_rng_fwrite(fp,r);
      fclose(fp);
      */
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
	      //std::cerr << t << '\t' << NSAM << '\t';
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
	      //std::cerr << X << endl;
	      ConditionalTraj(uni01,&path,N,s,dtp,ifreq);
	      //cerr << "sel phase..." << path.size() << "...";
	      selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,
			      &nlinks,&t,rho,nsites,N,path,X,1./double(2*N),1./double(4*k*N));
	      //cerr << "done" << endl;
	      neutral = true;
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,
					  sample_history,theta);
      PolySIM ad(&d);
      double simpi = (ad.ThetaPi()/double(nsites));
      if( fabs( obspi - simpi ) <= abs_tolerance ) ++prob;
      if( (1.-rel_tolerance)*obspi <= simpi && simpi <= (1.+rel_tolerance)*obspi ) ++prob2;
    }
  ofs.open(ofn,ios::app);
  ofs << 2*double(N)*s << '\t' << lambda << '\t' << theta << '\t' << rho << '\t'  
       << double(prob)/double(nreps) << '\t'
       << double(prob2)/double(nreps) << '\n';
  ofs.close();
}
