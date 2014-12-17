#include <Sequence/Coalescent/Coalescent.hpp>
#include <Sequence/PolySIM.hpp>
#include <sphase_twopop_nomig.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

enum EVENT {COAL1,COAL2,REC,SWEEP};

struct sortpos: public binary_function<polymorphicSite,polymorphicSite,bool>
{
  inline bool operator()(const polymorphicSite & lhs,
			 const polymorphicSite & rhs) const
  {
    return lhs.first <= rhs.first;
  }
};

int main( int argc, char **argv )
{
  if(argc != 12)
    {
      cerr << "rsweep_stochCG: recurrent hitch-hiking within the sampled region, and at linked sites (following Jensen et al. 2007, Genetics), with trajectories determined by the method of Coop and Griffiths\n"
	   << "usage: rsweep_linked_stoch nsam1 nsam2 nreps N s lambda rho_locus nsites theta_locus Tdiv seed\n"
	   << "where:\n"
	   << "\tlambda = expected # of sweeps per site per 4N generations\n"
	   << "\tnsam1 = sample size of neutral equilibrium pop\n"
	   << "\tnsam2 = sample size of population experiencing RHH\n"
	   << "\tTdiv = divergence time, in units of 4N generations\n";
      exit(1);
    }
  int argn=1;
  const int nsam1 = atoi(argv[argn++]);
  const int nsam2 = atoi(argv[argn++]);
  const int nreps = atoi(argv[argn++]);
  const int N = atoi(argv[argn++]);
  const double s = atof(argv[argn++]);
  const double lambda = atof(argv[argn++]);
  const double rho_locus = atof(argv[argn++]);
  const int nsites = atoi(argv[argn++]);
  const double theta = atof(argv[argn++]);
  const double Tdiv = atof(argv[argn++]);
  const unsigned seed = atoi(argv[argn++]);
  const double rho = rho_locus/double(nsites-1);

  if(lambda <= 0.) 
    {
      cerr << "error: lambda must be > 0!\n";
      exit(1);
    }
  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';
  vector<int> popconfig(1,nsam1);
  popconfig.push_back(nsam2);
  vector<chromosome> initialized_sample = init_sample( popconfig,
						       nsites );
  marginal initialized_marginal = init_marginal(nsam1+nsam2);
  vector<double> path;
  path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };

  double rcoal1,rcoal2,rrec,rsweepin=SEQMAXDOUBLE,rsweepout=SEQMAXDOUBLE,tmin,tcoal2,trec,t,tsweep;
  double maxd = 4.*double(N)*s;
  const int k = 50;
  const double dtp = 1./double(2*k*N),ifreq=1./double(2*N);
  int config[2];
  std::cerr << "//\n";
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      vector<double> fixedpos;
      bool merged = false;
      int NSAM = nsam1+nsam2;
      config[0]=nsam1;
      config[1]=nsam2;
      int nlinks = NSAM*(nsites-1);
      vector<chromosome> sample(initialized_sample);
      ARG sample_history(1,initialized_marginal);
      bool neutral = true;
      t=0.;
      while ( NSAM > 1 )
	{
	  //cerr << t << ' ' << config[0] << ' ' << config[1] << '\n';
	  if( neutral )
	    {
	      rcoal1 = (config[0]>1) ? double(config[0]*(config[0]-1)) : 0.;
	      rcoal2 = (config[1]>1) ? double(config[1]*(config[1]-1)) : 0.;
	      rrec = rho*double(nlinks);

	      //rate of sweeps within locus
	      rsweepin = lambda*double(nsites);
	      //rate of sweeps at linked sites
	      rsweepout = 2.*lambda*maxd/rho;

	      EVENT e = COAL1;
	      tmin = (rcoal1>0.) ? expo(1./rcoal1) : SEQMAXDOUBLE;
	      tcoal2 = (rcoal2>0.) ? expo(1./rcoal2) : SEQMAXDOUBLE;
	      if(tcoal2 <= tmin)
		{
		  tmin=tcoal2;
		  e = COAL2;
		}
	      trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	      if(trec <= tmin)
		{
		  tmin = trec;
		  e = REC;
		}
	      tsweep = (rsweepin+rsweepout>0. && config[1]>0) ? 
		expo(1./(rsweepin+rsweepout)): SEQMAXDOUBLE;
	      //cerr <<"tseep = " << tsweep << ", tmin = " << tmin << ", t = " << t << '\n';
	      if( tsweep <= tmin )
		{
		  //t += tsweep;
		  //neutral = false;
		  tmin = tsweep;
		  e = SWEEP;
		}

	      //cerr << e << ' ';
	      /*
	      if ( ((t + tmin >= Tdiv)
		    || (t < Tdiv && (config[0]==1&&config[1]==1)))
		   && ! merged )
	      */
	      if ( t+tmin >= Tdiv && !merged )
		{
		  t = Tdiv;
		  for(unsigned c=0;c<NSAM;++c)
		    {
		      sample[c].pop=1;
		    }
		  config[0]=0;
		  config[1]=NSAM;
		  merged = true;
		}
	      else if (e == SWEEP)
		{
		  t += tmin;
		  neutral = false;
		}
	      else // neutral phase events
		{
		  t += tmin;
		  if( e == COAL1 )
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							config[0],0);
		      int rv = coalesce(t,nsam1+nsam2,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      config[0] -= rv;
		      NSAM -= rv;
		    }
		  else if ( e == COAL2 )
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							config[1],1);
		      int rv = coalesce(t,nsam1+nsam2,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      config[1] -= rv;
		      NSAM -= rv;
		    }
		  else if ( e == REC )
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(uni01(),
								     nlinks,
								     sample.begin(),NSAM);
		      int pop = sample[pos_rec.first].pop;
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		      if(pop==0) config[0]++;
		      else if(pop==1) config[1]++;
		    }
// 		  if ( tcoal < trec ) //coalescence
// 		    {
// 		      std::pair<int,int> two = pick2(uni,NSAM);
// 		      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
// 				       &nlinks,&sample,&sample_history);
// 		    }
// 		  else //recombination
// 		    {
// 		      std::pair<int,int> pos_rec = pick_uniform_spot(gsl_rng_uniform(r),
// 								     nlinks,
// 								     sample.begin(),NSAM);
// 		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
// 					  &sample,&sample_history);
// 		      NSAM++;
// 		    }
		  if (NSAM < int(sample.size())/5)
		    {
		      sample.erase(sample.begin()+NSAM+1,sample.end());
		    }
		}
	    }
	  else //selective phase
	    {
	      std::cerr << t << '\t' ;//<< config[0] << ' ' << config[1] << '\t';
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
	      double xpos = double(X)/double(nsites);
	      if(xpos >= 0 && xpos < 1)
		{
		  fixedpos.push_back(xpos);
		}
	      std::cerr << X << '\n';
	      ConditionalTraj(std::bind(gsl_rng_uniform,r),&path,N,s,dtp,ifreq);
	      double told = t;
	      //cerr << "entering selective phase at time "<<t<<"...";
	      selective_phase_twopop_nomig(uni,uni01,sample,sample_history,nsam1+nsam2,&NSAM,
					   &nlinks,&t,&merged,rho,nsites,N,Tdiv,s,path,
					   X,k);//,1./double(2*N),1./double(4*k*N));
	      //cerr << "done"<<endl;
	      //cerr << t-told << '\n';
	      config[0]=config[1]=0;
	      for(unsigned c=0;c<NSAM;++c)
		{
		  config[sample[c].pop]++;
		}
	      neutral = true;
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,
					  sample_history,theta);
      if (!fixedpos.empty())
	{
	  vector<polymorphicSite> vps(d.sbegin(),d.send());
	  for(unsigned i=0;i<fixedpos.size();++i)
	    {
	      vps.push_back( make_pair(fixedpos[i],
				       string(nsam1,'0')+string(nsam2,'1')) );
	    }
	  sort(vps.begin(),vps.end(),sortpos());
	  d.assign(vps.begin(),vps.end());
	}

      cout << d << endl;
    }
}





