#pragma once

#include <string>
#include <cstdio>
#include <cstdlib>

namespace MLP {

    // Reproduces C#'s default double formatting: the shortest decimal
    // representation that round-trips back to the same double.
    inline std::string DblToStr(double v) {
        char buf[64];
        for (int p = 1; p <= 17; p++) {
            std::snprintf(buf, sizeof(buf), "%.*G", p, v);
            if (std::strtod(buf, nullptr) == v)
                break;
        }
        std::string s(buf);
        // C# writes exponents as E+16 / E-05, same as %G's E+16 / E-05 on glibc.
        return s;
    }

}
