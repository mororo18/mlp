#ifndef _ASSERT_HPP
#define _ASSERT_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace assertion {

    inline void error(const char *filename, int line) {
        std::cout << " Assertion failed. File " << filename << ":" << line << std::endl;
        std::exit(0);
    }

    inline void assert_eq(double a, double b, const char *filename, int line) {
        if (std::fabs(a - b) > std::numeric_limits<double>::epsilon()) {
            error(filename, line);
        }
    }

    inline void assert_eq(float a, float b, const char *filename, int line) {
        if (std::fabs(a - b) > std::numeric_limits<float>::epsilon()) {
            error(filename, line);
        }
    }

    inline void assert_eq(int a, int b, const char *filename, int line) {
        if (a != b) {
            error(filename, line);
        }
    }

    inline void assert_eq(bool a, bool b, const char *filename, int line) {
        if (a != b) {
            error(filename, line);
        }
    }

}

#endif
