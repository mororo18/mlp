#ifndef _TYPES_HPP
#define _TYPES_HPP

#include "preproc.hpp"

#include <cstddef>
#include <vector>

/*
 * Every array below keeps the original Fortran 1-based indexing:
 * the accessors translate the logical (1-based) subscripts into the
 * underlying 0-based storage, and the two dimensional ones keep the
 * column major layout used by Fortran.
 */

struct tMatrix {
    int n = 0;
    std::vector<typeReal> data;

    void allocate(int dimen) {
        n = dimen;
        data.assign(static_cast<std::size_t>(dimen) * static_cast<std::size_t>(dimen), typeReal(0.0));
    }

    inline typeReal &operator()(int i, int j) {
        return data[static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i - 1)];
    }

    inline const typeReal &operator()(int i, int j) const {
        return data[static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(n) + static_cast<std::size_t>(i - 1)];
    }
};

struct tInfo {
    int dimen = 0;
    tMatrix cost;

    static const int REINSERTION = 1;
    static const int OR_OPT_2    = 2;
    static const int OR_OPT_3    = 3;
    static const int SWAP        = 4;
    static const int TWO_OPT     = 5;

    typeReal fmax = static_cast<typeReal>(3.4028235E+38f);

    int reinsert_call = 0;

    std::vector<int> rnd;
    int rnd_index = 1;

    inline int rnd_at(int i) const { return rnd[static_cast<std::size_t>(i - 1)]; }
};

struct tSeqInfo {
    typeReal T, C, W;
};

struct tSolution {
    typeReal cost = typeReal(0.0);
    int s_size = 0;
    std::vector<int> s;
    std::vector<tSeqInfo> seq;
    int seq_dim = 0;

    inline int &s_at(int i) { return s[static_cast<std::size_t>(i - 1)]; }
    inline const int &s_at(int i) const { return s[static_cast<std::size_t>(i - 1)]; }

    inline tSeqInfo &seq_at(int j, int i) {
        return seq[static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(seq_dim) + static_cast<std::size_t>(j - 1)];
    }

    inline const tSeqInfo &seq_at(int j, int i) const {
        return seq[static_cast<std::size_t>(i - 1) * static_cast<std::size_t>(seq_dim) + static_cast<std::size_t>(j - 1)];
    }
};

#endif
