#ifndef _PREPROC
#define _PREPROC

#ifdef DOUBLE
#define typeReal 8
#else
#define typeReal 4
#endif

#define ASSERT(X,Y) call assert(X, Y, __FILE__, __LINE__)

#endif
