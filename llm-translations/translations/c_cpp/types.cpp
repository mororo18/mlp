#include "types.hpp"

#include <algorithm>

tSolution::tSolution(const tData & data) {
    size = sizeof(tSolution);
    s.assign(data.dimen + 1, 0);

    s_size = data.dimen + 1;

#ifdef MATRIX
    seq.assign(data.dimen + 1, std::vector<tSeq>(data.dimen + 1));

    size += s_size * sizeof(std::vector<tSeq>);
    size += (size_t)s_size * s_size * sizeof(tSeq);
#elif defined(FLAT)
    seq.assign((size_t)(data.dimen + 1) * (data.dimen + 1), tSeq());

    size += (size_t)s_size * s_size * sizeof(tSeq);
#endif

    MBsize = size / (1024.0 * 1024);

    cost = DBL_MAX;
}

void tSolution::cpy_to(tSolution & tgt, const tData & data) const {
    std::copy(s.begin(), s.begin() + (data.dimen + 1), tgt.s.begin());
    tgt.cost = cost;
}
