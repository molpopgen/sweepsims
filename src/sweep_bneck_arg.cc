#include <sweep_bneck_arg.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <Sequence/SeqConstants.hpp>
#include <util.hpp>
#include <sphase_bneck.hpp>
#include <iostream>
#include <cmath>
#include <functional>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;
using namespace Sequence;

//enum NEUT_EVENT {CA1,CA2,REC};

Sequence::arg sweep_bneck_arg(gsl_rng * r,
		    const int & n1,
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
		    const marginal & initialized_marginal,
		    const double & dt)
{
  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };

  const int nsam = n1;
  int NSAM = nsam;
  int nlinks = NSAM*(nsites-1);
  vector<chromosome> sample(initialized_sample);
  Sequence::arg sample_history(1,initialized_marginal);
  bool neutral = true;
  double t=0.;
  //double rc1,rc2,rrec,tc1,tc2,trec,tmin;
  double rcoal,rrec,tcoal,trec;
  while ( NSAM > 1 )
    {

      if( neutral )
	{
	  //cerr << t << ' ' << size << '\n';
	  rcoal = double(NSAM*(NSAM-1));
	  rrec = rho*double(nlinks);
	  tcoal = (t<tr||t>=tr+d) ? expo(1./rcoal) : expo(f/rcoal);
	  trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	  double tmin = min(tcoal,trec);
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
		  happens=false;
		}
	      else if ( t<(tr+d) && t+tmin >= (tr+d) )
		{
		  t=(tr+d);
		  happens=false;
		}
	      if(happens)
		{
		  t+=tmin;
		  if(tcoal < trec)
		    {
		      pair<int,int> two = pick2(uni,NSAM);
		      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		    }
		  else
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(uni01(),
								     nlinks,
								     sample.begin(),NSAM);
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		    }
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
	  const double told = t;
	  selective_phase(uni,uni01,sample,sample_history,
			  nsam,&NSAM,&nlinks,&t,f,rho,nsites,
			  N,s,path,X,dt);
	  //std::cerr << (t-told) << '\n';
	  
	  std::cerr <<(t-told)<< ' '
		    << -f*std::log(1./(2.*f*double(N)))/(2.*f*s*double(N)) << '\n';
	  
	  neutral = true;
	}
    }
  return sample_history;
}

Sequence::arg sweep_bneck_argCG(gsl_rng * r,
		      const int & n1,
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
		      const marginal & initialized_marginal,
		      const int & k,
		      const double & dt)
{
  std::function<double(void)> uni01 = [r](){ return gsl_rng_uniform(r); };
  std::function<double(const double&,const double&)> uni = [r](const double & a, const double & b){ return gsl_ran_flat(r,a,b); };
  std::function<double(const double&)> expo = [r](const double & mean){ return gsl_ran_exponential(r,mean); };
  std::function<double(const double&)> poiss = [r](const double & mean){ return gsl_ran_poisson(r,mean); };
  const int nsam = n1;
  int NSAM = nsam;
  int nlinks = NSAM*(nsites-1);
  vector<chromosome> sample(initialized_sample);
  Sequence::arg sample_history(1,initialized_marginal);
  bool neutral = true;
  double t=0.;
  //double rc1,rc2,rrec,tc1,tc2,trec,tmin;
  double rcoal,rrec,tcoal,trec;
  while ( NSAM > 1 )
    {

      if( neutral )
	{
	  //cerr << t << ' ' << size << '\n';
	  rcoal = double(NSAM*(NSAM-1));
	  rrec = rho*double(nlinks);
	  tcoal = (t<tr||t>=tr+d) ? expo(1./rcoal) : expo(f/rcoal);
	  trec = (rrec>0.) ? expo(1./rrec) : SEQMAXDOUBLE;
	  double tmin = min(tcoal,trec);
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
		  happens=false;
		}
	      else if ( t<(tr+d) && t+tmin >= (tr+d) )
		{
		  t=(tr+d);
		  happens=false;
		}
	      if(happens)
		{
		  t+=tmin;
		  if(tcoal < trec)
		    {
		      pair<int,int> two = pick2(uni,NSAM);
		      NSAM -= coalesce(t,nsam,NSAM,two.first,two.second,nsites,
					&nlinks,&sample,&sample_history);
		    }
		  else
		    {
		      std::pair<int,int> pos_rec = pick_uniform_spot(uni01(),
								     nlinks,
								     sample.begin(),NSAM);
		      nlinks -= crossover(NSAM,pos_rec.first,pos_rec.second,
					  &sample,&sample_history);
		      NSAM++;
		    }
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
	  double told=t;
	  selective_phaseCG(uni,uni01,sample,sample_history,
			    nsam,&NSAM,&nlinks,&t,f,rho,nsites,
			    N,s,path,X,k,dt);
	 
	  //	  std::cerr <<(t-told)<< '\n';
	  std::cerr <<(t-told)<< ' '
		    << -f*std::log(1./(2.*f*double(N)))/(2.*f*s*double(N)) << '\n';
	  neutral = true;
	}
    }
  return sample_history;
}

