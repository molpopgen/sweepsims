#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <sphase.hpp>
#include <util.hpp>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace std;
using namespace Sequence;

enum SEL_EVENT { FAVCOAL,UNFAVCOAL,FAVREC,FAVREC2UNFAV,
		 UNFAVREC,UNFAVREC2FAV,CHAOS };

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
  int pop = (event == FAVREC || event == FAVREC2UNFAV) ? 0 : 1;
  pair<int,int> two = pick_uniform_spot2(uni01(),nlinks_sel[pop],X,
					 sample.begin(),
					 *NSAM,pop);
  //if(sample[two.first].pop != pop) abort();
  int links0 = slinks(sample[two.first],X);
  int b = sample[two.first].first(),e=sample[two.first].last();
  if( two.second < b || two.second >= e )
    //recombination b/w ancestral material and selected site...
    {
      //if event is a "mover", then the sample config is changed,
      //else no effect
      if( event == FAVREC2UNFAV )
	{
	  sample[two.first].pop = 1;
	  config[0]--;
	  config[1]++;
	  nlinks_sel[0] -= links0;
	  nlinks_sel[1] += links0;
	}
      else if (event == UNFAVREC2FAV)
	{
	  sample[two.first].pop = 0;
	  config[0]++;
	  config[1]--;
	  nlinks_sel[0] += links0;
	  nlinks_sel[1] -= links0;
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
	  nlinks_sel[0] -= links0;
	  nlinks_sel[0] += slinks( sample[two.first], X );
	  nlinks_sel[0] += slinks( sample[(*NSAM)-1], X );
	}
      else if (event == UNFAVREC)
	{
	  nlinks_sel[1] -= links0;
	  nlinks_sel[1] += slinks( sample[two.first], X );
	  nlinks_sel[1] += slinks( sample[(*NSAM)-1], X );
	}
      else if (event == FAVREC2UNFAV)
	{
	  int migrant = (isleft) ? (*NSAM)-1 : two.first,
	    nonmigrant = (isleft) ? two.first :(*NSAM)-1 ;
	  sample[migrant].pop = 1;
	  config[0]--;
	  config[1]++;
	  nlinks_sel[0] -= links0;
	  nlinks_sel[0] += slinks( sample[nonmigrant], X );
	  nlinks_sel[1] += slinks( sample[migrant], X );
	}
      else if (event == UNFAVREC2FAV)
	{
	  int migrant = (isleft) ? two.first : (*NSAM)-1,
	    nonmigrant = (isleft) ? (*NSAM)-1 : two.first;
	  sample[migrant].pop = 0;
	  config[0]++;
	  config[1]--;
	  nlinks_sel[1] -= links0;
	  nlinks_sel[0] += slinks(sample[migrant], X );
	  nlinks_sel[1] += slinks(sample[nonmigrant], X );
	}
    }
#ifndef NDEBUG
  int t1=0,t2=0;
  for(int c=0 ; c < (*NSAM) ; ++c)
    {
      if(sample[c].pop==0) t1 += slinks(sample[c],X);
      else if (sample[c].pop==1) t2 += slinks(sample[c],X);
    }
  assert ( nlinks_sel[0] == t1 );
  assert ( nlinks_sel[1] == t2 );
#endif
}

void selective_phase( std::function<double(const double &, const double &)> & uni,
		      std::function<double()> & uni01,
		      vector<chromosome> & sample,
		      arg & sample_history,
		      const int & ttl_nsam,
		      int * NSAM,
		      int * nlinks,
		      double * t,
		      const double & f,
		      const double & rho,
		      const int & nsites,
		      const int & N,
		      const double & s,
		      const vector<double> & path,
		      const int & X,
		      const double & dt)
{
  //const double dt = 1./(2.*floor(f*double(N))),
  const double  pmin = 1./(2.*floor(f*double(N)));
  //const double tend = *t + f*(-std::log(pmin)/(2.*double(N)*f*s));

  int config[2], nlinks_sel[2];
  config[0] = *NSAM;
  config[1] = 0;
  nlinks_sel[0]=nlinks_sel[1]=0;
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
      nlinks_sel[0] += slinks(sample[c],X);
    }
  //fill arrays
  vector<double> pevents(6,0.);
  double tsel = 0.;
  int selgen = 0;
  //  while(config[0]>1)
  while(path[selgen]>=pmin)
    {
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      while( (r=uni01())== 1. );
      //selgen = int( tsel*4*N*f );
      //rejection method of braverman et al.
      while( (pr>r) && (path[selgen]>=pmin) )
	{
	  nevents=0;
	  for(unsigned pp = 0;pp<pevents.size();++pp) pevents[pp]=0.;
	  pevents[FAVCOAL] = (config[0]>1) ? dt*double(config[0]*(config[0]-1))/path[selgen] : 0.;
	  pevents[UNFAVCOAL] = (config[1]>1) ? dt*double(config[1]*(config[1]-1))/(1.-path[selgen]):0.;
	  pevents[FAVREC] = (rho>0.&&config[0]) ? dt*path[selgen]*rho*double(nlinks_sel[0]) : 0.;
	  pevents[FAVREC2UNFAV] = (rho>0.&&config[0]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[0]) : 0.;
	  pevents[UNFAVREC] = (rho>0.&&config[1]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[1]) : 0.;
	  pevents[UNFAVREC2FAV] = (rho>0.&&config[1]) ? dt*(path[selgen])*rho*double(nlinks_sel[1]) : 0.;
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
	  double ttemp = double(jump)*dt;
	  (*t) += f*ttemp;
	  tsel += ttemp;
	  SEL_EVENT event = CHAOS;
	  //figure out which event occurs
	  double rdm = uni01(),cum=0.;
	  //figure out which event occurs based on the 
	  //relative probability of each of the 6 types
	  transform(pevents.begin(),pevents.end(),pevents.begin(),
		    bind2nd(divides<double>(),sum));
	  for(int e = 0 ; e < 6 ; ++e)
	    {
	      cum += pevents[e];
	      if(rdm <= cum)
		{
		  event = SEL_EVENT(e);
		  e=6;
		}
	    }
	  assert(event != CHAOS);
	  
	  if( path[selgen] < pmin )
	    {
	      while( config[0] > 1 )
		{
		  pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,config[0],0);
		  int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				    nlinks,&sample,&sample_history);
		  (*NSAM) -= rv;
		  config[0]-=rv;
		}
	    }
	  else
	    {
	      if( event == FAVCOAL || event == UNFAVCOAL )
		{	  
		  int pop = (event == FAVCOAL) ? 0 : 1;
		  assert( *NSAM == config[0]+config[1] );
		  pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,
						    config[pop],pop);
		  int slinks_lost = slinks(sample[two.first],X);
		  int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				    nlinks,&sample,&sample_history);
		  (*NSAM) -= rv;
		  config[pop]-=rv;
		  assert( *NSAM == config[0]+config[1] );
		  bool flag = 0;
		  if (rv == 2)
		    {
		      flag=1;
		    }
		  else
		    {
		      slinks_lost -= slinks(sample[two.first],X);
		    }
		  slinks_lost += slinks(sample[*NSAM+rv-1-flag],X);
		  nlinks_sel[pop] -= slinks_lost;
#ifndef NDEBUG
		  //update nlinks_sel the hard way
		  //(note: figure out less dumb way to do this...)
		  
		  //nlinks_sel[0]=nlinks_sel[1]=0;
		  int t1=0,t2=0;
		  for(int c=0;c<(*NSAM);++c)
		    {
		      if(sample[c].pop==0)
			t1 +=  slinks(sample[c],X);
			//nlinks_sel[0] +=  slinks(sample[c],X);
		      else
			t2 += slinks(sample[c],X);
		      //nlinks_sel[1] +=  slinks(sample[c],X);
		    }
		  assert(t1 == nlinks_sel[0]);
		  assert(t2 == nlinks_sel[1]);
#endif
		}
	      else //some manner of recombination during the sweep.
		{
		  sel_rec(uni01,sample,sample_history,NSAM,
			  nlinks,nlinks_sel,config,event,X );
		}
	    }
	  assert (*NSAM == config[0]+config[1]);
	  if( (*NSAM) < int(sample.size())/5)
	    {
	      sample.erase(sample.begin()+(*NSAM)+1,sample.end());
	    }
#ifndef NDEBUG
	  int cn1=0,cn2=0;
	  for(int c = 0 ; c < (*NSAM) ;++c)
	    {
	      if (sample[c].pop==0) cn1++;
	      else if (sample[c].pop==1) cn2++;
	      else abort();
	    }
	  assert(cn1==config[0]);
	  assert(cn2==config[1]);
#endif
	}
    }
  //before exiting selective phase, put everyone back into same pop
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
    }
}

void selective_phaseCG( std::function<double(const double &, const double &)> & uni,
			std::function<double()> & uni01,
			vector<chromosome> & sample,
			arg & sample_history,
			const int & ttl_nsam,
			int * NSAM,
			int * nlinks,
			double * t,
			const double & f,
			const double & rho,
			const int & nsites,
			const int & N,
			const double & s,
			const vector<double> & path,
			const int & X,
			const int & k,
			const double & dt)
{
  //const double dt = 1./(2.*f*double(k*N)),
  const double pmin = 1./(2.*f*double(N));
  //const double tend = *t + -std::log(pmin)/(2.*double(N)*s);

  int config[2], nlinks_sel[2];
  config[0] = *NSAM;
  config[1] = 0;
  nlinks_sel[0]=nlinks_sel[1]=0;
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
      nlinks_sel[0] += slinks(sample[c],X);
    }
  //fill arrays
  vector<double> pevents(6,0.);
  double tsel = 0.;
  int selgen = 0;
  //  while(config[0]>1)
  //while(path[selgen]>=pmin)
  while(selgen < int(path.size()))
    {
      //std::cerr << selgen << '\n';
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      while( (r=uni01())== 1. );
      //selgen = int( tsel*4*N*k*f );
      //rejection method of braverman et al.
      assert(unsigned(selgen) < path.size());
      while( (pr>r) && (selgen<int(path.size())))//(path[selgen]>=pmin) )
	{
	  nevents=0;
	  for(unsigned pp = 0;pp<pevents.size();++pp) pevents[pp]=0.;
	  pevents[FAVCOAL] = (config[0]>1) ? dt*double(config[0]*(config[0]-1))/path[selgen] : 0.;
	  pevents[UNFAVCOAL] = (config[1]>1) ? dt*double(config[1]*(config[1]-1))/(1.-path[selgen]):0.;
	  pevents[FAVREC] = (rho>0.&&config[0]) ? dt*path[selgen]*rho*double(nlinks_sel[0]) : 0.;
	  pevents[FAVREC2UNFAV] = (rho>0.&&config[0]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[0]) : 0.;
	  pevents[UNFAVREC] = (rho>0.&&config[1]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[1]) : 0.;
	  pevents[UNFAVREC2FAV] = (rho>0.&&config[1]) ? dt*(path[selgen])*rho*double(nlinks_sel[1]) : 0.;
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
	  double ttemp = (double(jump)*dt);///2.;
	    (*t) += f*ttemp;
	    tsel += ttemp;
	    //if (! tsel <= double(path.size()-1)/double(4*N) )
	    //assert ( tsel <= double(path.size()-1)/double(4*N) );
	    SEL_EVENT event = CHAOS;
	    //figure out which event occurs
	    double rdm = uni01(),cum=0.;
	    //figure out which event occurs based on the 
	    //relative probability of each of the 6 types
	    transform(pevents.begin(),pevents.end(),pevents.begin(),
		      bind2nd(divides<double>(),sum));
	    for(int e = 0 ; e < 6 ; ++e)
	      {
		cum += pevents[e];
		if(rdm <= cum)
		  {
		    event = SEL_EVENT(e);
		    e=6;
		  }
	      }
	    assert(event != CHAOS);
	    if( event == FAVCOAL || event == UNFAVCOAL )
	      {	  
		int pop = (event == FAVCOAL) ? 0 : 1;
		assert( *NSAM == config[0]+config[1] );
		pair<int,int> two = pick2_in_deme(uni,sample,*NSAM,
						  config[pop],pop);
		//the use of slinks lost will make no sense
		//unless you study how Sequence::Coalesce works.
		//We need to keep track of the number of links
		//(relative to position of selected site) that
		//are lost from the sample due to the coalescent event.
		//We figure this out by doing what Sequence::Coalesce
		//does, but in reverse (ugh).
		int slinks_lost = slinks(sample[two.first],X);
		int rv = coalesce(*t,ttl_nsam,*NSAM,two.first,two.second,nsites,
				  nlinks,&sample,&sample_history);
		(*NSAM) -= rv;
		config[pop]-=rv;	
		assert( *NSAM == config[0]+config[1] );
		bool flag = 0;
		if (rv == 2)
		  {
		    flag=1;
		  }
		else
		  {
		    slinks_lost -= slinks(sample[two.first],X);
		  }
		slinks_lost += slinks(sample[*NSAM+rv-1-flag],X);
		nlinks_sel[pop] -= slinks_lost;

#ifndef NDEBUG
		//update nlinks_sel the hard way
		//(note: figure out less dumb way to do this...)
		//nlinks_sel[0]=nlinks_sel[1]=0;
		int t1=0,t2=0;
		for(int c=0;c<(*NSAM);++c)
		  {
		    if(sample[c].pop==0)
		      t1 +=  slinks(sample[c],X);
		    else
		      t2 +=  slinks(sample[c],X);
		  }
		assert( t1 == nlinks_sel[0] );
		assert( t2 == nlinks_sel[1] );
#endif
	      }
	    else //some manner of recombination during the sweep.
	      {
		sel_rec(uni01,sample,sample_history,NSAM,
			nlinks,nlinks_sel,config,event,X );
	      }
	}
      assert (*NSAM == config[0]+config[1]);
      if( (*NSAM) < int(sample.size())/5)
	{
	  sample.erase(sample.begin()+(*NSAM)+1,sample.end());
	}
#ifndef NDEBUG
      int cn1=0,cn2=0;
      for(int c = 0 ; c < (*NSAM) ;++c)
	{
	  if (sample[c].pop==0) cn1++;
	  else if (sample[c].pop==1) cn2++;
	  else assert(false);
	}
      assert(cn1==config[0]);
      assert(cn2==config[1]);
#endif
    }
  //}
  //before exiting selective phase, put everyone back into same pop
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
    }
}
