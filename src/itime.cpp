/*
* $URL$
* $Date$
* $Rev$
*/

#include "itime.hpp"

namespace LDPC
{

#ifdef WIN32

LARGE_INTEGER ITime::freq;
bool ITime::freqInitialized = QueryPerformanceFrequency(&ITime::freq);

#endif

};
