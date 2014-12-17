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
		 UNFAVREC,UNFAVREC2FAV,FAVREC2UNFAVMIG,UNFAVREC2FAVMIG,CHAOS };

void sel_rec( std::function<double(void)> & uni01, 
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

void selective_phase(std::function<double(const double &,const double &)> & uni,
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
		     //const double & s,
		     const vector<double> & path,
		     const int & X,
		     const double & pmin,
		     const double & dt,
		     const int & k)
{
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
  while(path[selgen]>=pmin)
    {
      //cerr << path[selgen] << '\n';
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      while( (r=uni01())== 1. );
      //rejection method of braverman et al.
      //if(selgen>=path.size())
      //cerr << selgen << ' ' << path.size() << ' ' << path[selgen] << '\n';
      //cerr << *NSAM << ' ' << config[0] << ' ' << config[1] << ' ' << selgen << ' ' << path.size() << ' ' 
      //<< path[selgen] << ' ' << pmin<< '\n';
      while( (pr>r) && (path[selgen]>=pmin) )
	{
	  //cerr << selgen << '\n';
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
      //cerr << "done\n";
      if(sum>0.)
	{
	  double ttemp = double(jump)*(dt);
	  (*t) += ttemp;
	  tsel += ttemp;
	  SEL_EVENT event = CHAOS;
	  //figure out which event occurs
	  double rdm = uni01(),cum=0.;
	  //cerr << "rdm = " << rdm << '\n';
	  //figure out which event occurs based on the 
	  //relative probability of each of the 6 types
	  transform(pevents.begin(),pevents.end(),pevents.begin(),
		    bind2nd(divides<double>(),sum));
	  for(int e = 0 ; e < 6 ; ++e)
	    {
	      //cerr << pevents[e] << ' ';
	      cum += pevents[e];
	      if(rdm <= cum)
		{
		  event = SEL_EVENT(e);
		  if(cum>0.)
		    //this is a guard against the RNG generating a 0,
		    //and us choosing an event with probability 0
		    {
		      break;
		    }
		}
	    }
	  assert(event != CHAOS);
	  //cerr << " -> "<< event << ' ' << (event == CHAOS) << '\n';
	  if( path[selgen] < pmin )//|| selgen >=  path.size())
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
		  //the use of slinks lost will make no sense
		  //unless you study how Coalesce works.
		  //We need to keep track of the number of links
		  //(relative to position of selected site) that
		  //are lost from the sample due to the coalescent event.
		  //We figure this out by doing what Coalesce
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
		  if( (*NSAM) < int(sample.size())/2)
		    {
		      sample.erase(sample.begin()+(*NSAM)+1,sample.end());
		    }
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
    }
  //before exiting selective phase, put everyone back into same pop
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
    }
  //cerr << selgen << "..." << path[selgen] << "..." << path[path.size()-1] << "...";
}

void selective_phase_partial( std::function<double(const double &,const double &)> & uni,
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
			      //const double & s,
			      const vector<double> & path,
			      const int & X,
			      const double & starting_freq,
			      const double & pmin,
			      const double & dt,
			      const int & k)
{
  int config[2], nlinks_sel[2];
  config[0] = 0;
  config[1] = 0;
  nlinks_sel[0]=nlinks_sel[1]=0;
  for(int c=0;c<(*NSAM);++c)
    {
      if(uni01()<=starting_freq)
	{
	  sample[c].pop=0;
	  nlinks_sel[0] += slinks(sample[c],X);
	  config[0]++;
	}
      else
	{
	  sample[c].pop=1;
	  nlinks_sel[1] += slinks(sample[c],X);
	  config[1]++;
	}
    }
  //fill arrays
  vector<double> pevents(6,0.);
  double tsel = 0.;
  int selgen = 0;
  while(path[selgen]>=pmin)
    {
      //cerr << path[selgen] << '\n';
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      while( (r=uni01())== 1. );
      //rejection method of braverman et al.
      //if(selgen>=path.size())
      //cerr << selgen << ' ' << path.size() << ' ' << path[selgen] << '\n';
      //cerr << *NSAM << ' ' << config[0] << ' ' << config[1] << ' ' << selgen << ' ' << path.size() << ' ' 
      //<< path[selgen] << ' ' << pmin<< '\n';
      while( (pr>r) && (path[selgen]>=pmin) )
	{
	  //cerr << selgen << '\n';
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
      //cerr << "done\n";
      if(sum>0.)
	{
	  double ttemp = double(jump)*(dt);
	  (*t) += ttemp;
	  tsel += ttemp;
	  SEL_EVENT event = CHAOS;
	  //figure out which event occurs
	  double rdm = uni01(),cum=0.;
	  //cerr << "rdm = " << rdm << '\n';
	  //figure out which event occurs based on the 
	  //relative probability of each of the 6 types
	  transform(pevents.begin(),pevents.end(),pevents.begin(),
		    bind2nd(divides<double>(),sum));
	  for(int e = 0 ; e < 6 ; ++e)
	    {
	      //cerr << pevents[e] << ' ';
	      cum += pevents[e];
	      if(rdm <= cum)
		{
		  event = SEL_EVENT(e);
		  if(cum>0.)
		    //this is a guard against the RNG generating a 0,
		    //and us choosing an event with probability 0
		    {
		      break;
		    }
		}
	    }
	  assert(event != CHAOS);
	  //cerr << " -> "<< event << ' ' << (event == CHAOS) << '\n';
	  if( path[selgen] < pmin )//|| selgen >=  path.size())
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
		  //the use of slinks lost will make no sense
		  //unless you study how Coalesce works.
		  //We need to keep track of the number of links
		  //(relative to position of selected site) that
		  //are lost from the sample due to the coalescent event.
		  //We figure this out by doing what Coalesce
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
		  if( (*NSAM) < int(sample.size())/2)
		    {
		      sample.erase(sample.begin()+(*NSAM)+1,sample.end());
		    }
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
    }
  //before exiting selective phase, put everyone back into same pop
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
    }
  //cerr << selgen << "..." << path[selgen] << "..." << path[path.size()-1] << "...";
}

void selective_phase_competing( std::function<double(const double &,const double &)> & uni,
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
				//const double & s,
				const vector<double> & path,
				const int & X,
				const double & pmin,
				const double & dt,
				const int & k)
{
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
  while(path[selgen]>=pmin && *NSAM > 1)
    {
      double pr=1.;
      double r,sum=0.;
      int nevents=0,jump=0;
      selgen = min(int(tsel*4*double(k*N)),int(path.size()-1));
      //selgen = int(tsel*2.*double(k*N));
      //selgen = min(int(2.*tsel*double(N)/double(k)),int(path.size()-1));
#ifndef NDEBUG
      /*
      if(selgen>=path.size())
	{
	  cerr <<"selgen = "<< selgen << ' ' << path.size() << ' ' << path[selgen] << ' ' 
	       << config[0] << ' ' << config[1] << ' ' << dt*double(config[0]*(config[0]-1))/path[selgen] << '\n';
	}
      */
      assert(selgen < path.size());
#endif
      //while( (r=uni01())== 1. );
      //rejection method of braverman et al.
      //      while( (pr>r) && (path[selgen]>=pmin) )
      //{
      //nevents=0;
      for(unsigned pp = 0;pp<pevents.size();++pp) pevents[pp]=0.;
      pevents[FAVCOAL] = (config[0]>1) ? dt*double(config[0]*(config[0]-1))/path[selgen] : 0.;
      pevents[UNFAVCOAL] = (config[1]>1) ? dt*double(config[1]*(config[1]-1))/(1.-path[selgen]):0.;
      pevents[FAVREC] = (rho>0.&&config[0]) ? dt*path[selgen]*rho*double(nlinks_sel[0]) : 0.;
      pevents[FAVREC2UNFAV] = (rho>0.&&config[0]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[0]) : 0.;
      pevents[UNFAVREC] = (rho>0.&&config[1]) ? dt*(1.-path[selgen])*rho*double(nlinks_sel[1]) : 0.;
      pevents[UNFAVREC2FAV] = (rho>0.&&config[1]) ? dt*(path[selgen])*rho*double(nlinks_sel[1]) : 0.;
	  /*
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
	  */
	  double ttemp = expo(1./pevents[0]);
	  SEL_EVENT event = SEL_EVENT(0);
	  for(int e = 1 ; e < 6 ; ++e)
	    {
	      double tm = expo(1./pevents[e]);
	      //cerr <<pevents[e] << ':'<< tm << ' ';
	      //if ( (!finite(ttemp) && finite(tm)) ||
	      //(finite(tm) && tm < ttemp) )
		if(tm < ttemp)
		{
		  ttemp=tm;
		  event = SEL_EVENT(e);
		}
	    }
	  (*t) += ttemp*dt;
	  tsel += ttemp*dt;
	  //cerr <<"\n:times become "<< *t << '\t' << tsel << '\n';
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
		  //the use of slinks lost will make no sense
		  //unless you study how Coalesce works.
		  //We need to keep track of the number of links
		  //(relative to position of selected site) that
		  //are lost from the sample due to the coalescent event.
		  //We figure this out by doing what Coalesce
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
		  if( (*NSAM) < int(sample.size())/2)
		    {
		      sample.erase(sample.begin()+(*NSAM)+1,sample.end());
		    }
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
  //before exiting selective phase, put everyone back into same pop
  for(int c=0;c<(*NSAM);++c)
    {
      sample[c].pop=0;
    }
}
