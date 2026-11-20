#ifndef _PREPROC_HPP
#define _PREPROC_HPP

#ifdef DOUBLE
typedef double typeReal;
#else
typedef float typeReal;
#endif

#define ASSERT(X, Y) assertion::assert_eq(X, Y, __FILE__, __LINE__)

#endif
