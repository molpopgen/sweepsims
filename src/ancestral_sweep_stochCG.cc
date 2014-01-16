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

enum EVENT { CA1, CA2, REC };

int main( int argc, char **argv )
{
  if( argc != 15 )
    {
      cerr << "ancestral_sweep_stochCG:\n"
	   << "usage: sweep nsam_ancestral nsam_derived nreps N tr tb f tau 2Ns X rho_site nsites theta_site seed\n"
	   << "where:\n"
	   << "\ttau is in units of 4N generations\n"
	   << "\tX is the position of the selected site, from 1 to nsites, inclusive\n";
      exit(1);
    }
  int argn=1;
  const int nsam = atoi(argv[argn++]);
  const int nsam2 = atoi(argv[argn++]);
  const int nreps = atoi(argv[argn++]);
  const int N = atoi(argv[argn++]);
  const double tr = atof(argv[argn++]);
  const double tb = atof(argv[argn++]);
  const double f = atof(argv[argn++]);
  const double tau = atof(argv[argn++]);
  const double alpha = atof(argv[argn++]);
  const int X = atoi(argv[argn++])-1; //user should input 1 to nsites
  const double rho = atof(argv[argn++]);
  const int nsites = atoi(argv[argn++]);
  const double theta = atof(argv[argn++]);
  const double s = alpha/(2.*N);
  const unsigned seed = atoi(argv[argn++]);
  //cerr << rho << ' ' << nsites << ' ' << theta << '\n';
  if ( tb <= tr )
    {
      cerr << "error: tb must be > tr\n";
      exit(1);
    }
  if(tau<=tb)
    {
      cerr << "error: tau must be > tb\n";
      exit(1);
    }
  copy(argv,argv+argc,ostream_iterator<char *>(cout," "));
  cout << '\n';
  //const double rho = irho*double(nsites)/double(nsites-1);
  vector<int> config(1,nsam);
  config.push_back(nsam2);
  vector<chromosome> initialized_sample = init_sample( config,
						       nsites );
  marginal initialized_marginal = init_marginal(nsam+nsam2);

  vector<double> path;
  path.reserve(int(std::floor(4*std::log(2.*double(N))/s)));

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  //double ar;
  gsl_uniform01 uni01(r); 
  gsl_uniform uni(r);     
  gsl_exponential expo(r);
  gsl_poisson poiss(r); 

  //double rcoal,rrec,tcoal,trec,t;
  const int k = 10;
  const double dtp = 1./double(2*k*N),ifreq=1./double(2*N);
  double t,tmin,ttemp,fcurr=1.,sum,den,rdm;
  EVENT e;
  double rates[3];
  for(int rep = 0 ; rep < nreps ; ++rep)
    {
      int NSAM = nsam+nsam2;
      int N1 = nsam, N2 = nsam2;
      int nlinks = NSAM*(nsites-1);
      vector<chromosome> sample(initialized_sample);
      Sequence::arg sample_history(1,initialized_marginal);
      bool neutral = true;
      t=0.;
      bool swept=false,merged=false;
      while ( NSAM > 1 )
	{
	  //cerr << t << ' '<< swept << ' ' << neutral << ' ' << NSAM << ' ' << N1 << ' ' << N2 << '\n';
	  /*
	  if(swept)
	    {
	      cerr << N1*(N1-1) << ' ' << rho*double(nlinks) << '\n';
	    }
	  */
#ifndef NDEBUG
	  int check1=0,check2=0;
	  for(int i=0;i<NSAM;++i)
	    {
	      if( sample[i].pop == 0 ) ++check1;
	      else ++check2;
	    }
	  assert( check1==N1 );
	  assert( check2==N2 );
	  assert(N1+N2==NSAM);
#endif
	  if( neutral )
	    {
	      fcurr = ( t >= tr && t < tb ) ? f : 1.; 
	      rates[CA1] = (N1>1) ? double(N1*(N1-1)) : 0.;
	      rates[CA2] = (N2>1) ? double(N2*(N2-1))/fcurr : 0.;
	      rates[REC] = (rho>0&&nlinks>0) ? rho*double(nlinks) : 0.;
	      den=rates[CA1]+rates[CA2]+rates[REC];
	      sum=0.;
	      rdm = gsl_rng_uniform(r);
	      for(unsigned i=0;i<3;++i)
		{
		  sum += rates[i];
		  if ( rdm <= sum/den )
		    {
		      e = EVENT(i);
		      i=3;
		    }
		}
	      tmin = expo(1./rates[unsigned(e)]);

	      assert(isfinite(tmin));

	      bool event = true;

	      if(!merged && (t<tb&&t+tmin>=tb))
		{
		  assert(!swept);
		  t=tb;
		  for(int i=0;i<NSAM;++i)
		    {
		      sample[i].pop=0;
		    }
		  N1+=N2;
		  N2=0;
		  assert(N1==NSAM);
		  merged = true;
		  event = false;
		}
	      if(event && !swept && t<tau && t+tmin>tau)
		{
		  t=tau;
		  event  = false;
		  neutral=false;
		}		
#ifndef NDEBUG
	      int checkl=0;
	      for(int i=0;i<NSAM;++i) checkl += sample[i].links();
	      assert(checkl==nlinks);
#endif
	      if(event)
		{
		  t += tmin;
		  if ( e==CA1 ) //coalescence
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							N1,0);
		      assert(two.first != two.second);
		      assert( sample[two.first].pop == 0);
		      assert( sample[two.second].pop == 0);
		      int rv = coalesce(t,nsam+nsam2,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      N1 -= rv;
		      NSAM -= rv;
		      assert(N1+N2==NSAM);
		    }
		  else if (e == CA2)
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							N2,1);
		      assert(two.first != two.second);
		      assert( sample[two.first].pop == 1);
		      assert( sample[two.second].pop == 1);
		      int rv = coalesce(t,nsam+nsam2,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      N2 -= rv;
		      NSAM -= rv;
		      assert(N1+N2==NSAM);
		    }
		  else if(e==REC)
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(gsl_rng_uniform(r),
								     nlinks,
								     sample.begin(),NSAM);
		      int pop = sample[pos_rec.first].pop;
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		      if(pop==0) N1++;
		      else N2++;
		      assert(N1+N2==NSAM);
		    }
		}
	      if (NSAM < int(sample.size())/5)
		{
		  sample.erase(sample.begin()+NSAM+1,sample.end());
		}   
	    }
	  else if (!neutral)//selective phase
	    {
#ifndef NDEBUG
	      for(int i=0;i<NSAM;++i) 
		{
		  assert ( sample[i].first() >= 0 &&
			   sample[i].last() < nsites );
		  assert ( sample[i].pop == 0 );
		}
#endif
	      assert(!swept);
	      assert(t==tau);
	      assert(N2==0);
	      assert(N1+N2==NSAM);
	      int sel_site = (X>=0) ? X : int(gsl_ran_flat(r,0.,nsites));
	      cerr  << sel_site << '\n';
	      ConditionalTraj(uni01,&path,N,s,dtp,ifreq);
	      //cerr << t << ' ' << NSAM << ' ';

	      selective_phase(uni,uni01,sample,sample_history,nsam+nsam2,&NSAM,
			      &nlinks,&t,rho,nsites,N,path,X,
			      1./double(2*N),1./double(4*k*N),k);

	      //cerr << t << ' ' << NSAM << '\n';
	      N1=NSAM;
	      swept=true;
	      neutral = true;
#ifndef NDEBUG
	      for(int i=0;i<NSAM;++i) 
		{
		  assert ( sample[i].first() >= 0 &&
			   sample[i].last() < nsites );
		  assert ( sample[i].pop == 0 );
		}
#endif
	    }
	}
      minimize_arg(&sample_history);
      SimData d = infinite_sites_sim_data(poiss,uni,nsites,sample_history,theta*double(nsites));
      cout << d << '\n' ;
    }
}
