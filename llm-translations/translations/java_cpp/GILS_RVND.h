#ifndef GILS_RVND_H
#define GILS_RVND_H

#include "tArray.h"
#include "tInfo.h"
#include "tSolution.h"

#include <memory>

class GILS_RVND {
private:
    // parameters
    int                     Iils;
    static constexpr int    Imax = 10;
    static constexpr double R[] = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                                   0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    static constexpr int    R_size = 26;

    bool                    verbose;

public:
    std::unique_ptr<tInfo> info;

    GILS_RVND(bool verbose);

private:
    void subseq_load(tSolution & solut, tInfo & info);
    tArray construction(double alpha, tInfo & info);
    bool search_swap(tSolution & solut, tInfo & info);
    bool search_two_opt(tSolution & solut, tInfo & info);
    bool search_reinsertion(tSolution & solut, tInfo & info, int opt);
    void RVND(tSolution & solut, tInfo & info);
    tArray perturb(tArray & sl, tInfo & info);

public:
    void solve();
};

#endif
