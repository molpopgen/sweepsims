#include <util.hpp>

int slinks( const Sequence::chromosome & c,
		   const int & X )
{
  return ( std::max(X,c.last()) - std::min(X,c.first()) );
}

std::pair<int,int> 
pick_uniform_spot2(const double & random_01,
		   const int & nlinks,
		   const int & X,
		   std::vector<Sequence::chromosome>::const_iterator sample_begin,
		   const unsigned & current_nsam,
		   const int & deme)
{
  int pos = int(random_01*double(nlinks))+1;
  int recombinant = 0,len=0;
  while( (sample_begin+recombinant) <
	 (sample_begin+current_nsam) )
    {
      if( (sample_begin+recombinant)->pop == deme )
	{
	  len = slinks( *(sample_begin+recombinant),X );
	  if(pos <= len)break;
	  pos -= len;
	}
      recombinant++;
    }
  assert( (sample_begin+recombinant)->pop == deme );
  int rpos = std::min(X,(sample_begin+recombinant)->begin()->beg) + pos - 1;
  return std::make_pair(recombinant,rpos);
}

std::pair<int,int> 
pick_uniform_spot3(const double & random_01,
		   const int & nlinks,
		   std::vector<Sequence::chromosome>::const_iterator sample_begin,
		   const unsigned & current_nsam,
		   const int & deme)
{
  int pos = int(random_01*double(nlinks))+1;
  int recombinant = 0,len=0;
  while( (sample_begin+recombinant) <
	 (sample_begin+current_nsam) )
    {
      if( (sample_begin+recombinant)->pop == deme )
	{
	  //len = slinks( *(sample_begin+recombinant),X );
	  len = (sample_begin+recombinant)->links();
	  if(pos <= len)break;
	  pos -= len;
	}
      recombinant++;
    }
  assert( (sample_begin+recombinant)->pop == deme );
  //int rpos = std::min(X,(sample_begin+recombinant)->begin()->beg) + pos - 1;
  int rpos = (sample_begin+recombinant)->first()+pos-1;
  return std::make_pair(recombinant,rpos);
}
