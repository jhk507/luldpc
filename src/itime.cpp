/*
* $URL$
* $Date$
* $Rev$
*/

#include "itime.hpp"

#ifdef WIN32

LARGE_INTEGER ITime::freq;
bool ITime::freqInitialized = QueryPerformanceFrequency(&ITime::freq);

#endif
