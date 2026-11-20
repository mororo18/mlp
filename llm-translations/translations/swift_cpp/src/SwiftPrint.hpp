#pragma once

#include <string>
#include <vector>

// Reproduces Swift's default textual representation of a Double
// (shortest round-tripping form, always with a decimal point or exponent).
std::string describeDouble(double v);

// Reproduces Swift's default textual representation of an array of Int.
std::string describeSequence(const std::vector<int> &v);
