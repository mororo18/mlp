#pragma once

#include <cstddef>
#include <vector>

enum Operator {
    Swap = 0,
    Reinsertion = 1,
    OrOpt2 = 2,
    OrOpt3 = 3,
    TwoOpt = 4,
};

struct SubseqInfo {
    double t;
    double c;
    double w;

    static inline SubseqInfo zeored() { return SubseqInfo{0.0, 0.0, 0.0}; }
};

class SubseqMatrix {
#ifdef FLAT
    std::size_t row_size;
    std::vector<SubseqInfo> data;
#else
    std::vector<std::vector<SubseqInfo>> data;
#endif

  public:
    SubseqMatrix() = default;

    explicit SubseqMatrix(std::size_t n)
#ifdef FLAT
        : row_size(n + 1), data((n + 1) * (n + 1), SubseqInfo::zeored()) {
    }
#else
        : data(n + 1, std::vector<SubseqInfo>(n + 1, SubseqInfo::zeored())) {
    }
#endif

    inline const SubseqInfo &operator()(std::size_t row, std::size_t col) const {
#ifdef FLAT
        return data[row_size * row + col];
#else
        return data[row][col];
#endif
    }

    inline SubseqInfo &operator()(std::size_t row, std::size_t col) {
#ifdef FLAT
        return data[row_size * row + col];
#else
        return data[row][col];
#endif
    }
};

class CostMatrix {
#ifdef FLAT
    std::size_t row_size;
    std::vector<double> data;
#else
    std::vector<std::vector<double>> data;
#endif

  public:
    CostMatrix() = default;

    explicit CostMatrix(std::size_t n)
#ifdef FLAT
        : row_size(n), data(n * n, 0.0) {
    }
#else
        : data(n, std::vector<double>(n, 0.0)) {
    }
#endif

    inline const double &operator()(std::size_t row, std::size_t col) const {
#ifdef FLAT
        return data[row_size * row + col];
#else
        return data[row][col];
#endif
    }

    inline double &operator()(std::size_t row, std::size_t col) {
#ifdef FLAT
        return data[row_size * row + col];
#else
        return data[row][col];
#endif
    }
};

struct Data {
    CostMatrix c;
    std::size_t dimen;
    std::vector<std::size_t> rnd;
    std::size_t rnd_index;
};

struct Solution {
    SubseqMatrix seq;
    std::vector<std::size_t> s;
    double cost;

    inline std::size_t operator[](std::size_t index) const { return s[index]; }
};
