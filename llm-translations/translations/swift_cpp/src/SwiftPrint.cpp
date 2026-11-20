#include "SwiftPrint.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

std::string describeDouble(double v) {
	if (std::isnan(v)) {
		return "nan";
	}
	if (std::isinf(v)) {
		return v > 0 ? "inf" : "-inf";
	}
	if (v == 0.0) {
		return std::signbit(v) ? "-0.0" : "0.0";
	}

	char buf[64];
	for (int p = 1; p <= 17; p++) {
		std::snprintf(buf, sizeof(buf), "%.*e", p - 1, v);
		if (std::strtod(buf, nullptr) == v) {
			break;
		}
	}

	std::string s(buf);
	size_t epos = s.find('e');
	int exp10 = std::atoi(s.c_str() + epos + 1);
	std::string mant = s.substr(0, epos);

	bool neg = false;
	if (!mant.empty() && mant[0] == '-') {
		neg = true;
		mant.erase(0, 1);
	}

	std::string digits;
	for (char c : mant) {
		if (c != '.') {
			digits += c;
		}
	}
	while (digits.size() > 1 && digits.back() == '0') {
		digits.pop_back();
	}

	std::string out;
	if (exp10 >= -4 && exp10 < 17) {
		if (exp10 >= 0) {
			if ((int)digits.size() > exp10 + 1) {
				out = digits.substr(0, exp10 + 1) + "." + digits.substr(exp10 + 1);
			} else {
				out = digits + std::string(exp10 + 1 - digits.size(), '0') + ".0";
			}
		} else {
			out = "0." + std::string(-exp10 - 1, '0') + digits;
		}
	} else {
		out = digits.substr(0, 1);
		if (digits.size() > 1) {
			out += "." + digits.substr(1);
		}
		char eb[16];
		std::snprintf(eb, sizeof(eb), "e%c%02d", exp10 < 0 ? '-' : '+', std::abs(exp10));
		out += eb;
	}

	return (neg ? "-" : "") + out;
}

std::string describeSequence(const std::vector<int> &v) {
	std::string out = "[";
	for (size_t i = 0; i < v.size(); i++) {
		if (i != 0) {
			out += ", ";
		}
		out += std::to_string(v[i]);
	}
	out += "]";
	return out;
}
