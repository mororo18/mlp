#pragma once

#include "subseq.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace data {

struct LoadResult {
    std::size_t dimension;
    CostMatrix c;
    std::vector<std::size_t> rnd;
};

LoadResult load(const std::string &filename);

} // namespace data
