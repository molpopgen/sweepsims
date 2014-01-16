#include <Sequence/Coalescent/Trajectories.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/bind.hpp>
using namespace Sequence;
using namespace std;

int main(int argc, char **argv)
{
  const unsigned N = atoi(argv[1]);
  const double s = atof(argv[2]);
  const unsigned seed =atoi(argv[3]);

  gsl_rng * r =  gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  vector<double>path;
  ConditionalTraj(boost::bind(gsl_rng_uniform,r),&path,N,s,1./(4*N),1/double(2*N));
  copy(path.begin(),path.end(),ostream_iterator<double>(cout,"\n"));
}
