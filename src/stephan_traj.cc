#include <detpath.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>

using namespace std;

int main(int argc, char **argv)
{
  const unsigned N = atoi(argv[1]);
  const double s = atof(argv[2]);

  vector<double>path = detpath(N,s);
  copy(path.begin(),path.end(),ostream_iterator<double>(cout,"\n"));
}
