#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstring>
#include <cfloat>
#include <cstdlib>
#include <vector>

using Real = double;

struct tData {
    std::vector<std::vector<Real>> cost;
    int dimen = 0;
    std::vector<int> rnd;
    int rnd_index = 0;
};

enum tNeighborhood {
    REINSERTION = 1,
    OR_OPT_2 = 2,
    OR_OPT_3 = 3,
    SWAP = 4,
    TWO_OPT = 5
};

struct tSeqInfo {
    Real T = 0.0, C = 0.0, W = 0.0;
};

using tSeq = tSeqInfo;

struct tSolution {
#ifdef MATRIX
    std::vector<std::vector<tSeq>> seq;
#elif defined(FLAT)
    std::vector<tSeq> seq;
#endif
    //std::vector<std::vector<std::vector<Real>>> seq;
    Real cost = DBL_MAX;
    std::vector<int> s;
    int s_size = 0;
    size_t size = 0;
    float MBsize = 0.0f;

    tSolution() = default;
    explicit tSolution(const tData & data);

    void cpy_to(tSolution & tgt, const tData & data) const;

#ifdef FLAT
    int to_1D(const int i, const int j, const int size) const;
#endif

    /*==========================SET==========================*/

    void seq_set_C(int i, int j, Real value);
    void seq_set_T(int i, int j, Real value);
    void seq_set_W(int i, int j, Real value);

    /*==========================GET==========================*/

    Real seq_get_C(int i, int j) const;
    Real seq_get_T(int i, int j) const;
    Real seq_get_W(int i, int j) const;
};

/*==========================inline==========================*/

#ifdef FLAT
inline
int tSolution::to_1D(const int i, const int j, const int size) const {
    return i * size + j;
}
#endif


inline
void tSolution::seq_set_C(int i, int j, Real value) {
#ifdef MATRIX
    seq[i][j].C = value;
#elif defined(FLAT)
    seq[to_1D(i, j, s_size)].C = value;
#endif
}

inline
void tSolution::seq_set_T(int i, int j, Real value) {
#ifdef MATRIX
    seq[i][j].T = value;
#elif defined(FLAT)
    seq[to_1D(i, j, s_size)].T = value;
#endif
}

inline
void tSolution::seq_set_W(int i, int j, Real value) {
#ifdef MATRIX
    seq[i][j].W = value;
#elif defined(FLAT)
    seq[to_1D(i, j, s_size)].W = value;
#endif
}


inline
Real tSolution::seq_get_C(int i, int j) const {
#ifdef MATRIX
    return seq[i][j].C;
#elif defined(FLAT)
    return seq[to_1D(i, j, s_size)].C;
#endif
}

inline
Real tSolution::seq_get_T(int i, int j) const {
#ifdef MATRIX
    return seq[i][j].T;
#elif defined(FLAT)
    return seq[to_1D(i, j, s_size)].T;
#endif
}

inline
Real tSolution::seq_get_W(int i, int j) const {
#ifdef MATRIX
    return seq[i][j].W;
#elif defined(FLAT)
    return seq[to_1D(i, j, s_size)].W;
#endif
}

#endif
