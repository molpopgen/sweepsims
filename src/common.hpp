#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <config.h>

#include <Sequence/Coalescent/Coalescent.hpp>

#ifdef HAVE_SEQUENCE_COALSIM
using ARG = Sequence::coalsim::arg;
using namespace Sequence;
using namespace Sequence::coalsim;
#else
using ARG = Sequence::arg;
using namespace Sequence;
#endif

#endif
