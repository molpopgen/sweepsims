#include <sweep_der_arg.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <Sequence/SeqConstants.hpp>
#include <util.hpp>
#include <sphase_der.hpp>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

enum NEUT_EVENT {CA1,CA2,REC};

Sequence::arg sweep_der_arg(gsl_rng * r,
		  const int & n1,
		  const int & n2,
		  const int & nsites,
		  const double & rho,
		  const double & tr,
		  const double & d,
		  const double & f,
		  const int & N,
		  const double & s,
		  const double & tau,
		  const int & X,
		  const vector<double> & path,
		  const vector<chromosome> & initialized_sample,
		  const marginal & initialized_marginal)
{
  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };
  const int nsam = n1+n2;
  int NSAM = nsam,N1=n1,N2=n2;
  int nlinks = NSAM*(nsites-1);
  vector<chromosome> sample(initialized_sample);
  Sequence::arg sample_history(1,initialized_marginal);
  bool neutral = true;
  double t=0.;
  double size = 1.;
  double rc1,rc2,rrec,tc1,tc2,trec,tmin;
  while ( NSAM > 1 )
    {
      //cerr << t << ' ' << neutral << ' ' << N1 << ' ' << N2 << '\n';
      if( neutral )
	{
	  //cerr << t << ' ' << size << '\n';
	  rc1 = (N1>0) ? double(N1*(N1-1)) : 0.;
	  rc2 = (N2>0) ? double(N2*(N2-1)) : 0.;
	  rrec = rho*double(nlinks);
	  
	  tc1 = (rc1>0) ? ((t>=tr+d) ? expo( (1.+f)/rc1 ) : expo(1./rc1)) : SEQMAXDOUBLE;
	  tmin = tc1;
	  NEUT_EVENT event = CA1;
	  tc2 = expo(size/rc2);
	  if(tc2<tmin)
	    {
	      tmin=tc2;
	      event = CA2;
	    }
	  trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	  if (trec < tmin)
	    {
	      tmin = trec;
	      event = REC;
	    }
	  if ( (t==0.&&tau==0.) || (t < tau && t+tmin >= tau) )
	    {
	      t = tau;
	      neutral = false;
	    }
	  else // neutral phase events
	    {
	      bool happens = true;
	      
	      if( t < tr && t+tmin >= tr )
		{
		  t=tr;
		  size=f;
		  happens=false;
		}
	      else if ( t<(tr+d) && t+tmin >= (tr+d) )
		{
		  t=(tr+d);
		  size=1;
		  happens=false;
		  //merge all pops into pop 0
		  for(int c=0;c<NSAM;++c)
		    {
		      sample[c].pop=0;
		    }
		  N1 += N2;
		  N2 = 0;
		  assert(NSAM==N1);
		}
	      if(happens)
		{
		  t += tmin;
		  if (event == CA1)//coal. in ancestrap pop
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							N1,0);
		      int rv = coalesce(t,nsam,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      N1 -= rv;
		      NSAM -= rv;
		    }
		  else if (event == CA2) //coal in derived pop
		    {
		      pair<int,int> two = pick2_in_deme(uni,sample,NSAM,
							N2,1);
		      int rv = coalesce(t,nsam,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		      N2 -= rv;
		      NSAM -= rv;
		    }
		  else if (event == REC) //recombination
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(uni01(),
								     nlinks,
								     sample.begin(),NSAM);
		      int pop = sample[pos_rec.first].pop;
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		      if(pop==0) N1++;
		      else if(pop==1) N2++;
		    }
		  assert(NSAM==N1+N2);
		}
	      if (NSAM < int(sample.size())/5)
		{
		  sample.erase(sample.begin()+NSAM+1,sample.end());
		}
	    }
	}
      else //selective phase
	{
	  //see sphase_der.hpp and sphase_der.cc
	  //	  cerr << "phase entered at " << t << '\n';
	  selective_phase(uni,uni01,sample,sample_history,nsam,&NSAM,
			  &nlinks,&t,rho,nsites,N,f,s,path,X);
	  //selective_phase_competing(uni,uni01,expo,sample,sample_history,nsam,&NSAM,
	  //&nlinks,&t,rho,nsites,N,f,s,path,X);
	  N1=N2=0;
	  for(int c=0;c<NSAM;++c)
	    {
	      if(sample[c].pop==0) N1++;
	      else if(sample[c].pop==1) N2++;
	    }
	  neutral = true;
	  if( t > tr+d ) 
	    {
	      cerr << " t = "<< t << '\n';
	      abort();
	    }
	}
    }
  return sample_history;
}
