#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <sphase.hpp>
#include <util.hpp>
#include <functional>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Sequence;

//these are the 3 demes we need to deal with
const int ANCESTRAL = 0; //ancestral pop
const int FAVORED = 1;   //favored class w/in derived pop
const int UNFAVORED = 2; //unfavored class w/in derived pop

//the types of events that can happen
enum SEL_EVENT { FAVCOAL,UNFAVCOAL,FAVREC,FAVREC2UNFAV,
		 UNFAVREC,UNFAVREC2FAV,ANCCOAL,ANCREC,CHAOS };

//in case we want to output
ostream & operator<<( ostream & o, const SEL_EVENT & e )
{
  if ( e == FAVCOAL )
    o << "FAVCOAL";
  else if ( e == UNFAVCOAL )
    o << "UNFAVCOAL";
  else if ( e == FAVREC )
    o << "FAVREC";
  else if ( e == UNFAVREC )
    o << "UNFAVREC";
  else if ( e == FAVREC2UNFAV )
    o << "FAVREC2UNFAV";
  else if ( e == UNFAVREC2FAV )
    o << "UNFAVREC2FAV";
  else if ( e == ANCCOAL )
    o << "ANCCOAL";
  else if ( e == ANCREC )
    o << "ANCREC";
  return o;
}

void sel_rec( std::function<double()> & uni01, 
	      vector<chromosome> & sample,
	      arg & sample_history,
	      int * NSAM,
	      int * nlinks,
	      int nlinks_sel[],
	      int config[],
	      const SEL_EVENT & event,
	      const int & X)
{
  int pop = (event == FAVREC || event == FAVREC2UNFAV) ? FAVORED : UNFAVORED;
  pair<int,int> two = pick_uniform_spot2(uni01(),nlinks_sel[pop],X,
					 sample.begin(),
					 *NSAM,pop);

  int links0 = slinks(sample[two.first],X);
  int b = sample[two.first].first(),e=sample[two.first].last();
  if( two.second < b || two.second >= e )
    //recombination b/w ancestral material and selected site...
    {
      //if event is a "mover", then the sample config is changed,
      //else no effect
      if( event == FAVREC2UNFAV )
	{
	  sample[two.first].pop = UNFAVORED;
	  config[FAVORED]--;
	  config[UNFAVORED]++;
	  nlinks_sel[FAVORED] -= links0;
	  nlinks_sel[UNFAVORED] += links0;
	}
      else if (event == UNFAVREC2FAV)
	{
	  sample[two.first].pop = FAVORED;
	  config[FAVORED]++;
	  config[UNFAVORED]--;
	  nlinks_sel[FAVORED] += links0;
	  nlinks_sel[UNFAVORED] -= links0;
	}
    }
  else
    {
      (*nlinks) -= crossover(*NSAM,two.first,two.second,
			     &sample,&sample_history);
      (*NSAM)++;
      config[pop]++;
      //is selected site to the left of the position of xover?
      bool isleft = (X <= two.second);
      
      if (event == FAVREC)
	{
	  nlinks_sel[FAVORED] -= links0;
	  nlinks_sel[FAVORED] += slinks( sample[two.first], X );
	  nlinks_sel[FAVORED] += slinks( sample[(*NSAM)-1], X );
	}
      else if (event == UNFAVREC)
	{
	  nlinks_sel[UNFAVORED] -= links0;
	  nlinks_sel[UNFAVORED] += slinks( sample[two.first], X );
	  nlinks_sel[UNFAVORED] += slinks( sample[(*NSAM)-1], X );
	}
      else if (event == FAVREC2UNFAV)
	{
	  int migrant = (isleft) ? (*NSAM)-1 : two.first,
	    nonmigrant = (isleft) ? two.first :(*NSAM)-1 ;
	  sample[migrant].pop = UNFAVORED;
	  config[FAVORED]--;
	  config[UNFAVORED]++;
	  nlinks_sel[FAVORED] -= links0;
	  nlinks_sel[FAVORED] += slinks( sample[nonmigrant], X );
	  nlinks_sel[UNFAVORED] += slinks( sample[migrant], X );
	}
      else if (event == UNFAVREC2FAV)
	{
	  int migrant = (isleft) ? two.first : (*NSAM)-1,
	    nonmigrant = (isleft) ? (*NSAM)-1 : two.first;
	  sample[migrant].pop = FAVORED;
	  config[FAVORED]++;
	  config[UNFAVORED]--;
	  nlinks_sel[UNFAVORED] -= links0;
	  nlinks_sel[FAVORED] += slinks(sample[migrant], X );
	  nlinks_sel[UNFAVORED] += slinks(sample[nonmigrant], X );
	}
    }
#ifndef NDEBUG
  int t0=0,t1=0,t2=0;
  for(int c=0 ; c < (*NSAM) ; ++c)
    {
      if(sample[c].pop==ANCESTRAL) t0 += sample[c].links();
      else if(sample[c].pop==FAVORED) t1 += slinks(sample[c],X);
      else if (sample[c].pop==UNFAVORED) t2 += slinks(sample[c],X);
    }
  assert ( nlinks_sel[ANCESTRAL] == t0 );
  assert ( nlinks_sel[FAVORED] == t1 );
  assert ( nlinks_sel[UNFAVORED] == t2 );
#endif
}

void selective_phase( std::function<double(const double&,const double&)> & uni,
		      std::function<double()> & uni01,
		      vector<chromosome> & sample,
		      arg & sample_history,
		      const int & ttl_nsam,
		      int * NSAM,
		      int * nlinks,
		      double * t,
		      const double & rho,
		      const int & nsites,
		      const int & N,
		      const double & f,
		      const double & s,
		      const vector<double> & path,
		      const int & X)
{
  //const double oldt=*t;
  const double dt = 1./(4.*floor(f*double(N))), pmin = 1./(2.*floor(f*double(N))),dtN = 1./(4.*double(N));
  double tsel = 0.,accumt=0.;
  //const double tend = *t + -std::log(pmin)/(2.*double(f*N)*s);

  int config[3], nlinks_sel[3];
  config[ANCESTRAL]=config[FAVORED]=config[UNFAVORED]=
    nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
  for(int c=0;c<(*NSAM);++c)
    {
      if(sample[c].pop == 0)
	{
	  //chromo is in ancestral population
	  config[ANCESTRAL]++;
	  nlinks_sel[ANCESTRAL] += sample[c].links();
	}
      else if (sample[c].pop == 1)
	{
	  //sample in bottlenecked population
	  config[FAVORED]++;
	  nlinks_sel[FAVORED] += slinks(sample[c],X);
	}
    }
  assert( config[ANCESTRAL]+config[FAVORED] == *NSAM );
  assert( config[UNFAVORED] == 0 );
  //fill arrays
  vector<double> pevents(8,0.);
  int selgen = 0;
  while(path[selgen]>=pmin)
    {
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      while( (r=uni01())== 1. );
      selgen = int(tsel*4.*f*double(N));
      //cerr << selgen << '\n';
      //rejection method of braverman et al.
      while( (pr>r) && (path[selgen]>=pmin) )
	{
	  nevents=0;
	  for(unsigned pp = 0;pp<pevents.size();++pp) pevents[pp]=0.;
	  pevents[FAVCOAL] = (config[FAVORED]>1) ? dt*double(config[FAVORED]*(config[FAVORED]-1))/path[selgen] : 0.;
	  pevents[UNFAVCOAL] = (config[UNFAVORED]>1) ? dt*double(config[UNFAVORED]*(config[UNFAVORED]-1))/(1.-path[selgen]):0.;
	  pevents[FAVREC] = (rho>0.&&config[FAVORED]) ? dt*path[selgen]*rho*double(nlinks_sel[FAVORED]) : 0.;
	  pevents[FAVREC2UNFAV] = (rho>0.&&config[FAVORED]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[FAVORED]) : 0.;
	  pevents[UNFAVREC] = (rho>0.&&config[UNFAVORED]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[UNFAVORED]) : 0.;
	  pevents[UNFAVREC2FAV] = (rho>0.&&config[UNFAVORED]) ? dt*(path[selgen])*rho*double(nlinks_sel[UNFAVORED]) : 0.;
	  pevents[ANCCOAL] = (config[ANCESTRAL]>1) ? dtN*double(config[ANCESTRAL]*(config[ANCESTRAL]-1)) : 0.;
	  pevents[ANCREC] = (rho>0.&&config[ANCESTRAL]) ? dtN*rho*double(nlinks_sel[ANCESTRAL]) : 0.;

	  //hypotheticals
	  //compare the rate of coalescence in the favored class, to that in the ancestral class,
	  //assuming the sampple sizes are the same at the current size
	  //if we have done things right, the rate in the ancestral class will be f*path[selgen]
	  //lower than than in the favored class...
	  /*
	    cerr << dt*double(config[FAVORED]*(config[FAVORED]-1))/path[selgen] << ' '
	    << dtN*double(config[FAVORED]*(config[FAVORED]-1)) << ' ' << f << ' ' << path[selgen] << '\n';
	  */
	  sum = accumulate(pevents.begin(),pevents.end(),0.);
	  if(sum>=1.)
	    {
	      nevents++; //more than 1 event in a generation
	    }
	  
	  pr *= (1.-sum);
	  if(nevents==0) 
	    { 
	      ++selgen;
	      ++jump;
	    }
	}
      if(sum>0.)
	{
	  //Accounting for 2 different time scales here!
	  double ttemp = double(jump)*dtN;
	  //increment time on the sweep in units of 4N generations
	  (*t) += ttemp;  
	  accumt += ttemp;
	  //increment time along the trajectory of the beneficial allele
	  //in units of 4Nf generations.
	  tsel += double(jump)*dt;
	  SEL_EVENT event = CHAOS;
	  double rdm = uni01(),cum=0.;
	  //figure out which event occurs based on the 
	  //relative probability of each of the 8 types
	  transform(pevents.begin(),pevents.end(),pevents.begin(),
		    bind2nd(divides<double>(),sum));
	  for(int e = 0 ; e < 8 ; ++e)
	    {
	      cum += pevents[e];
	      if(rdm < cum)
		{
		  event = SEL_EVENT(e);
		  e=8;
		}
	    }
	  assert(event != CHAOS);
	  
	  /*
	  cerr << *t << ' ' << tsel << ' ' << config[FAVORED] << ' ' << pevents[FAVCOAL] << ' ' 
	       << config[ANCESTRAL] << ' ' << pevents[ANCCOAL] << ' ' << event << '\n';	  
	  */
	  if( path[selgen] < pmin )
	    {
	      while( config[FAVORED] > 1 )
		{
		  pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,config[FAVORED],FAVORED);
		  int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				    nlinks,&sample,&sample_history);
		  (*NSAM) -= rv;
		  config[FAVORED]-=rv;
		}
	      nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
	      for(int c=0;c<(*NSAM);++c)
		{
		  if(sample[c].pop == ANCESTRAL)
		    nlinks_sel[ANCESTRAL] += sample[c].links();
		  else if(sample[c].pop==FAVORED)
		    nlinks_sel[FAVORED] +=  slinks(sample[c],X);
		  else if (sample[c].pop == UNFAVORED )
		    nlinks_sel[UNFAVORED] +=  slinks(sample[c],X);
		}
	    }
	  else
	    {
	      if( event == FAVCOAL || event == UNFAVCOAL || event == ANCCOAL)
		{	  
		  int pop = FAVORED;
		  if(event == UNFAVCOAL) pop = UNFAVORED;
		  else if (event == ANCCOAL) pop = ANCESTRAL;
		  assert( *NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED] );
		  pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,
						    config[pop],pop);
		  assert(sample[two.first].pop == sample[two.second].pop);
		  assert(sample[two.first].pop == pop);
		  int slinks_lost = (pop!=ANCESTRAL) ? slinks(sample[two.first],X) : sample[two.first].links();
		  int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				    nlinks,&sample,&sample_history);
		  (*NSAM) -= rv;
		  config[pop]-=rv;
		  assert( *NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED] );
		  bool flag = 0;
		  if (rv == 2)
		    {
		      flag=1;
		    }
		  else
		    {
		      slinks_lost -= (pop!=ANCESTRAL)?slinks(sample[two.first],X):sample[two.first].links();
		    }
		  slinks_lost += (pop!=ANCESTRAL)?slinks(sample[*NSAM+rv-1-flag],X):sample[*NSAM+rv-1-flag].links();
		  nlinks_sel[pop] -= slinks_lost;

#ifndef NDEBUG
		  //update nlinks_sel the hard way
		  //(note: figure out less dumb way to do this...)
		  
		  //nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
		  int t1=0,t2=0,t3=0;
		  for(int c=0;c<(*NSAM);++c)
		    {
		      if(sample[c].pop == ANCESTRAL)
			t1+= sample[c].links();
			//nlinks_sel[ANCESTRAL] += sample[c].links();
		      else if(sample[c].pop==FAVORED)
			t2+= slinks(sample[c],X);
			//nlinks_sel[FAVORED] +=  slinks(sample[c],X);
		      else if (sample[c].pop == UNFAVORED )
			t3+= slinks(sample[c],X);
			//nlinks_sel[UNFAVORED] +=  slinks(sample[c],X);
		    }
		  assert( t1 == nlinks_sel[0] );
		  assert( t2 == nlinks_sel[1] );
		  assert( t3 == nlinks_sel[2] );
#endif
		}
	      else if (event == ANCREC)
		{
		  pair<int,int> two = pick_uniform_spot3(uni01(),nlinks_sel[ANCESTRAL],
							 sample.begin(),
							 *NSAM,ANCESTRAL);
		  assert(sample[two.first].pop==ANCESTRAL);
		  int rv = crossover(*NSAM,two.first,two.second,
				     &sample,&sample_history);
		  //update sample sizes
		  (*NSAM)++;
		  config[ANCESTRAL]++;
		  //update various links, etc.
		  (*nlinks) -= rv;
		  nlinks_sel[ANCESTRAL] -= rv;
		}
	      else //some manner of recombination in the derived pop during the sweep.
		{
		  sel_rec(uni01,sample,sample_history,NSAM,
			  nlinks,nlinks_sel,config,event,X );
		}
	    }
	  assert (*NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED]);
	  if( (*NSAM) < int(sample.size())/5)
	    {
	      sample.erase(sample.begin()+(*NSAM)+1,sample.end());
	    }
#ifndef NDEBUG
	  //check chromo configs
	  int cn0=0,cn1=0,cn2=0;
	  int nl=0,n0=0,n1=0,n2=0;
	  for(int c = 0 ; c < (*NSAM) ;++c)
	    {
	      nl += sample[c].links();
	      if(sample[c].pop == ANCESTRAL)
		{
		  cn0++;
		  n0 += sample[c].links();
		}
	      else if (sample[c].pop==FAVORED) 
		{
		  cn1++;
		  n1 += slinks(sample[c],X);
		}
	      else if (sample[c].pop==UNFAVORED) 
		{
		  cn2++;
		  n2 += slinks(sample[c],X);
		}
	      else abort();
	    }
	  assert(cn0==config[ANCESTRAL]);
	  assert(cn1==config[FAVORED]);
	  assert(cn2==config[UNFAVORED]);
	  assert(nl == *nlinks);
	  assert(n0 == nlinks_sel[ANCESTRAL]);
	  assert(n1 == nlinks_sel[FAVORED]);
	  assert(n2 == nlinks_sel[UNFAVORED]);
	  
	  //check nlinks stuff
	  
#endif
	}
    }
  //before exiting selective phase, put everyone back into appropriate pop
  for(int c=0;c<(*NSAM);++c)
    {
      if(sample[c].pop == ANCESTRAL)
	{
	  //no change
	}
      else //goes into bnecked population
	{
	  sample[c].pop = 1;
	  //cerr << sample[c].first() << ' ' << sample[c].last() << '\n';
	}
    }
  
  /*
  cerr << "tsel = "<< tsel << ' ' << *t << ' ' << 2.*f*double(N)*s << ' ' << -log(1./double(2*f*N))/(2*double(f*N)*s)<< ' ' 
       << f*(-log(1./double(2*floor(double(N)*f)))/(2*double(N)*f*s)) 
       << ' ' << selgen << ' ' << path.size() << '\n'; 
  cerr << "accumulated t/4N0 = " << accumt << '\n';
  cerr << "oldt = " << oldt << ", newt = " << *t << '\n';
  */
  //exit(1);
}

void selective_phase_competing( std::function<double(const double&,const double&)> & uni,
				std::function<double()> & uni01,
				std::function<double(const double &)> & expo,
				vector<chromosome> & sample,
				arg & sample_history,
				const int & ttl_nsam,
				int * NSAM,
				int * nlinks,
				double * t,
				const double & rho,
				const int & nsites,
				const int & N,
				const double & f,
				const double & s,
				const vector<double> & path,
				const int & X)
/*!
  Uses competing exponentials to determine times.  Much better for v. large regions
*/
{
  //const double oldt=*t;
  const double dt = 1./(4.*floor(f*double(N))), pmin = 1./(2.*floor(f*double(N))),dtN = 1./(4.*double(N));
  double tsel = 0.,accumt=0.;
  //const double tend = *t + -std::log(pmin)/(2.*double(f*N)*s);

  int config[3], nlinks_sel[3];
  config[ANCESTRAL]=config[FAVORED]=config[UNFAVORED]=
    nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
  for(int c=0;c<(*NSAM);++c)
    {
      if(sample[c].pop == 0)
	{
	  //chromo is in ancestral population
	  config[ANCESTRAL]++;
	  nlinks_sel[ANCESTRAL] += sample[c].links();
	}
      else if (sample[c].pop == 1)
	{
	  //sample in bottlenecked population
	  config[FAVORED]++;
	  nlinks_sel[FAVORED] += slinks(sample[c],X);
	}
    }
  assert( config[ANCESTRAL]+config[FAVORED] == *NSAM );
  assert( config[UNFAVORED] == 0 );
  //fill arrays
  vector<double> pevents(8,0.);
  int selgen = 0;

  while(path[selgen]>=pmin)
    {
      double r;//,sum=0.;
      while( (r=uni01())== 1. );
      selgen = int(tsel*4.*f*double(N));

      for(unsigned pp = 0;pp<pevents.size();++pp) pevents[pp]=0.;
      pevents[FAVCOAL] = (config[FAVORED]>1) ? dt*double(config[FAVORED]*(config[FAVORED]-1))/path[selgen] : 0.;
      pevents[UNFAVCOAL] = (config[UNFAVORED]>1) ? dt*double(config[UNFAVORED]*(config[UNFAVORED]-1))/(1.-path[selgen]):0.;
      pevents[FAVREC] = (rho>0.&&config[FAVORED]) ? dt*path[selgen]*rho*double(nlinks_sel[FAVORED]) : 0.;
      pevents[FAVREC2UNFAV] = (rho>0.&&config[FAVORED]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[FAVORED]) : 0.;
      pevents[UNFAVREC] = (rho>0.&&config[UNFAVORED]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[UNFAVORED]) : 0.;
      pevents[UNFAVREC2FAV] = (rho>0.&&config[UNFAVORED]) ? dt*(path[selgen])*rho*double(nlinks_sel[UNFAVORED]) : 0.;
      pevents[ANCCOAL] = (config[ANCESTRAL]>1) ? dtN*double(config[ANCESTRAL]*(config[ANCESTRAL]-1)) : 0.;
      pevents[ANCREC] = (rho>0.&&config[ANCESTRAL]) ? dtN*rho*double(nlinks_sel[ANCESTRAL]) : 0.;
      
      double ttemp = expo(1./pevents[0]);
      SEL_EVENT event = SEL_EVENT(0);
      for(int e = 1 ; e < 8 ; ++e)
	{
	  double tm = expo(1./pevents[e]);
	  if(tm < ttemp)
	    {
	      ttemp=tm;
	      event = SEL_EVENT(e);
	    }
	}
      
      //Accounting for 2 different time scales here!
      //increment time on the sweep in units of 4N generations
      (*t) += ttemp*dtN;  
      accumt += ttemp*dtN;
      //increment time along the trajectory of the beneficial allele
      //in units of 4Nf generations.
      tsel += ttemp*dt;//double(jump)*dt;
      
      if( path[selgen] < pmin )
	{
	  while( config[FAVORED] > 1 )
	    {
	      pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,config[FAVORED],FAVORED);
	      int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				nlinks,&sample,&sample_history);
	      (*NSAM) -= rv;
	      config[FAVORED]-=rv;
	    }
	  nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
	  for(int c=0;c<(*NSAM);++c)
	    {
	      if(sample[c].pop == ANCESTRAL)
		nlinks_sel[ANCESTRAL] += sample[c].links();
	      else if(sample[c].pop==FAVORED)
		nlinks_sel[FAVORED] +=  slinks(sample[c],X);
	      else if (sample[c].pop == UNFAVORED )
		nlinks_sel[UNFAVORED] +=  slinks(sample[c],X);
	    }
	}
      else
	{
	  if( event == FAVCOAL || event == UNFAVCOAL || event == ANCCOAL)
	    {	  
	      int pop = FAVORED;
	      if(event == UNFAVCOAL) pop = UNFAVORED;
	      else if (event == ANCCOAL) pop = ANCESTRAL;
	      assert( *NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED] );
	      pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,
						config[pop],pop);
	      assert(sample[two.first].pop == sample[two.second].pop);
	      assert(sample[two.first].pop == pop);
	      int slinks_lost = (pop!=ANCESTRAL) ? slinks(sample[two.first],X) : sample[two.first].links();
	      int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				nlinks,&sample,&sample_history);
	      (*NSAM) -= rv;
	      config[pop]-=rv;
	      assert( *NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED] );
	      bool flag = 0;
	      if (rv == 2)
		{
		  flag=1;
		}
	      else
		{
		  slinks_lost -= (pop!=ANCESTRAL)?slinks(sample[two.first],X):sample[two.first].links();
		}
	      slinks_lost += (pop!=ANCESTRAL)?slinks(sample[*NSAM+rv-1-flag],X):sample[*NSAM+rv-1-flag].links();
	      nlinks_sel[pop] -= slinks_lost;

#ifndef NDEBUG
	      //update nlinks_sel the hard way
	      //(note: figure out less dumb way to do this...)
		  
	      //nlinks_sel[ANCESTRAL]=nlinks_sel[FAVORED]=nlinks_sel[UNFAVORED]=0;
	      int t1=0,t2=0,t3=0;
	      for(int c=0;c<(*NSAM);++c)
		{
		  if(sample[c].pop == ANCESTRAL)
		    t1+= sample[c].links();
		  //nlinks_sel[ANCESTRAL] += sample[c].links();
		  else if(sample[c].pop==FAVORED)
		    t2+= slinks(sample[c],X);
		  //nlinks_sel[FAVORED] +=  slinks(sample[c],X);
		  else if (sample[c].pop == UNFAVORED )
		    t3+= slinks(sample[c],X);
		  //nlinks_sel[UNFAVORED] +=  slinks(sample[c],X);
		}
	      assert( t1 == nlinks_sel[0] );
	      assert( t2 == nlinks_sel[1] );
	      assert( t3 == nlinks_sel[2] );
#endif
	    }
	  else if (event == ANCREC)
	    {
	      pair<int,int> two = pick_uniform_spot3(uni01(),nlinks_sel[ANCESTRAL],
						     sample.begin(),
						     *NSAM,ANCESTRAL);
	      assert(sample[two.first].pop==ANCESTRAL);
	      int rv = crossover(*NSAM,two.first,two.second,
				 &sample,&sample_history);
	      //update sample sizes
	      (*NSAM)++;
	      config[ANCESTRAL]++;
	      //update various links, etc.
	      (*nlinks) -= rv;
	      nlinks_sel[ANCESTRAL] -= rv;
	    }
	  else //some manner of recombination in the derived pop during the sweep.
	    {
	      sel_rec(uni01,sample,sample_history,NSAM,
		      nlinks,nlinks_sel,config,event,X );
	    }
	}
      assert (*NSAM == config[ANCESTRAL]+config[FAVORED]+config[UNFAVORED]);
      if( (*NSAM) < int(sample.size())/5)
	{
	  sample.erase(sample.begin()+(*NSAM)+1,sample.end());
	}
#ifndef NDEBUG
      //check chromo configs
      int cn0=0,cn1=0,cn2=0;
      int nl=0,n0=0,n1=0,n2=0;
      for(int c = 0 ; c < (*NSAM) ;++c)
	{
	  nl += sample[c].links();
	  if(sample[c].pop == ANCESTRAL)
	    {
	      cn0++;
	      n0 += sample[c].links();
	    }
	  else if (sample[c].pop==FAVORED) 
	    {
	      cn1++;
	      n1 += slinks(sample[c],X);
	    }
	  else if (sample[c].pop==UNFAVORED) 
	    {
	      cn2++;
	      n2 += slinks(sample[c],X);
	    }
	  else abort();
	}
      assert(cn0==config[ANCESTRAL]);
      assert(cn1==config[FAVORED]);
      assert(cn2==config[UNFAVORED]);
      assert(nl == *nlinks);
      assert(n0 == nlinks_sel[ANCESTRAL]);
      assert(n1 == nlinks_sel[FAVORED]);
      assert(n2 == nlinks_sel[UNFAVORED]);
	  
      //check nlinks stuff
	  
#endif
    }
  //}
  //before exiting selective phase, put everyone back into appropriate pop
  for(int c=0;c<(*NSAM);++c)
    {
      if(sample[c].pop == ANCESTRAL)
	{
	  //no change
	}
      else //goes into bnecked population
	{
	  sample[c].pop = 1;
	  //cerr << sample[c].first() << ' ' << sample[c].last() << '\n';
	}
    }
  
  /*
    cerr << "tsel = "<< tsel << ' ' << *t << ' ' << 2.*f*double(N)*s << ' ' << -log(1./double(2*f*N))/(2*double(f*N)*s)<< ' ' 
    << f*(-log(1./double(2*floor(double(N)*f)))/(2*double(N)*f*s)) 
    << ' ' << selgen << ' ' << path.size() << '\n'; 
    cerr << "accumulated t/4N0 = " << accumt << '\n';
    cerr << "oldt = " << oldt << ", newt = " << *t << '\n';
  */
  //exit(1);
}
