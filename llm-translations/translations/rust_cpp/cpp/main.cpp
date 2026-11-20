#include "data.hpp"
#include "subseq.hpp"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <limits>
#include <string>
#include <vector>

static const double F64_EPSILON = std::numeric_limits<double>::epsilon();
static const double F64_MAX = std::numeric_limits<double>::max();

// Mimics Rust's `{}` (Display) formatting for f64: shortest representation
// that round-trips, always in decimal notation.
static std::string fmt_f64(double v) {
    char buf[1200];
    auto res = std::to_chars(buf, buf + sizeof(buf), v, std::chars_format::fixed);
    return std::string(buf, res.ptr);
}

// Mimics Rust's `{:?}` (Debug) formatting for f64: same as Display but always
// with a fractional part.
static std::string fmt_f64_debug(double v) {
    std::string s = fmt_f64(v);
    if (s.find('.') == std::string::npos) {
        s += ".0";
    }
    return s;
}

// Mimics Rust's `{:?}` formatting for Vec<usize>.
static std::string fmt_vec(const std::vector<std::size_t> &v) {
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); i++) {
        if (i != 0) {
            s += ", ";
        }
        s += std::to_string(v[i]);
    }
    s += "]";
    return s;
}

static void update_subseq_info_matrix(Solution &solut, const Data &info) {
    for (std::size_t i = 0; i <= info.dimen; i++) {
        const int k = 1 - static_cast<int>(i) - (i == 0 ? 1 : 0);

        solut.seq(i, i).t = 0.0;
        solut.seq(i, i).c = 0.0;
        solut.seq(i, i).w = (i != 0) ? 1.0 : 0.0;

        for (std::size_t j = i + 1; j <= info.dimen; j++) {
            const std::size_t j_prev = j - 1;

            const double t =
                info.c(solut[j_prev], solut[j]) + solut.seq(i, j_prev).t;
            solut.seq(i, j).t = t;

            const double c = solut.seq(i, j).t + solut.seq(i, j_prev).c;
            solut.seq(i, j).c = c;

            const double w = static_cast<double>(static_cast<int>(j) + k);
            solut.seq(i, j).w = w;
        }
    }

    solut.cost = solut.seq(0, info.dimen).c;
}

static long partition(std::vector<std::size_t> &arr, long left, long right,
                      const Data &info, std::size_t r) {
    const std::size_t pivot = arr[right];
    long i = left - 1;
    for (long j = left; j < right; j++) {
        if (info.c(r, arr[j]) < info.c(r, pivot)) {
            i += 1;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

static void quicksort(std::vector<std::size_t> &arr, long left, long right,
                      const Data &info, std::size_t r) {
    if (left < right) {
        const long pivot = partition(arr, left, right, info, r);
        quicksort(arr, left, pivot - 1, info, r);
        quicksort(arr, pivot + 1, right, info, r);
    }
}

static void sort(std::vector<std::size_t> &arr, std::size_t r,
                 const Data &info) {
    quicksort(arr, 0, static_cast<long>(arr.size()) - 1, info, r);
}

static std::vector<std::size_t> construction(double alpha, Data &info) {
    (void)alpha;
    std::vector<std::size_t> s(1, 0);

    std::vector<std::size_t> c_list;
    for (std::size_t i = 1; i < info.dimen; i++) {
        c_list.push_back(i);
    }

    std::size_t r = 0;
    while (!c_list.empty()) {
        sort(c_list, r, info);

        const std::size_t index = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        const std::size_t c = c_list[index];
        r = c;
        c_list.erase(c_list.begin() + index);
        s.push_back(c);
    }
    s.push_back(0);

    return s;
}

static void swap(std::vector<std::size_t> &s, std::size_t i, std::size_t j) {
    std::swap(s[i], s[j]);
}

static void reverse(std::vector<std::size_t> &s, std::size_t i, std::size_t j) {
    std::reverse(s.begin() + i, s.begin() + j + 1);
}

static void reinsert(std::vector<std::size_t> &s, std::size_t i, std::size_t j,
                     std::size_t pos) {
    if (i < pos) {
        for (std::size_t k = i; k <= j; k++) {
            const std::size_t tmp = s[i];
            s.insert(s.begin() + pos, tmp);
            s.erase(s.begin() + i);
        }
    } else {
        for (std::size_t k = i; k <= j; k++) {
            const std::size_t tmp = s[j];
            s.erase(s.begin() + j);
            s.insert(s.begin() + pos, tmp);
        }
    }
}

static bool search_swap(Solution &solut, const Data &info) {
    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    double cost_concat_4;

    double cost_best = F64_MAX;
    double cost_new;
    std::size_t I = 0;
    std::size_t J = 0;

    for (std::size_t i = 1; i < info.dimen - 1; i++) {
        const std::size_t i_prev = i - 1;
        const std::size_t i_next = i + 1;

        cost_concat_1 =
            solut.seq(0, i_prev).t + info.c(solut[i_prev], solut[i_next]);
        cost_concat_2 = cost_concat_1 + solut.seq(i, i_next).t +
                        info.c(solut[i], solut[i_next + 1]);

        cost_new = solut.seq(0, i_prev).c +
                   solut.seq(i, i_next).w * cost_concat_1 +
                   info.c(solut[i_next], solut[i]) +
                   solut.seq(i_next + 1, info.dimen).w * cost_concat_2 +
                   solut.seq(i_next + 1, info.dimen).c;

        if (cost_new < cost_best) {
            cost_best = cost_new - F64_EPSILON;
            I = i;
            J = i_next;
        }

        for (std::size_t j = i_next + 1; j < info.dimen; j++) {
            const std::size_t j_next = j + 1;
            const std::size_t j_prev = j - 1;

            cost_concat_1 =
                solut.seq(0, i_prev).t + info.c(solut[i_prev], solut[j]);
            cost_concat_2 = cost_concat_1 + info.c(solut[j], solut[i_next]);
            cost_concat_3 = cost_concat_2 + solut.seq(i_next, j_prev).t +
                            info.c(solut[j_prev], solut[i]);
            cost_concat_4 = cost_concat_3 + info.c(solut[i], solut[j_next]);

            cost_new = solut.seq(0, i_prev).c + cost_concat_1 +
                       solut.seq(i_next, j_prev).w * cost_concat_2 +
                       solut.seq(i_next, j_prev).c + cost_concat_3 +
                       solut.seq(j_next, info.dimen).w * cost_concat_4 +
                       solut.seq(j_next, info.dimen).c;

            if (cost_new < cost_best) {
                cost_best = cost_new - F64_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - F64_EPSILON) {
        swap(solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    } else {
        return false;
    }
}

static bool search_two_opt(Solution &solut, const Data &info) {
    double cost_new;
    double cost_best = F64_MAX;

    double cost_concat_1;
    double cost_concat_2;

    std::size_t I = 0;
    std::size_t J = 0;

    for (std::size_t i = 1; i < info.dimen - 1; i++) {
        const std::size_t i_prev = i - 1;
        double rev_seq_cost = solut.seq(i, i + 1).t;

        for (std::size_t j = i + 2; j < info.dimen; j++) {
            const std::size_t j_next = j + 1;

            rev_seq_cost += info.c(solut[j - 1], solut[j]) *
                            (solut.seq(i, j).w - 1.0);

            cost_concat_1 =
                solut.seq(0, i_prev).t + info.c(solut[j], solut[i_prev]);
            cost_concat_2 = cost_concat_1 + solut.seq(i, j).t +
                            info.c(solut[j_next], solut[i]);

            cost_new = solut.seq(0, i_prev).c +
                       solut.seq(i, j).w * cost_concat_1 + rev_seq_cost +
                       solut.seq(j_next, info.dimen).w * cost_concat_2 +
                       solut.seq(j_next, info.dimen).c;

            if (cost_new < cost_best) {
                cost_best = cost_new - F64_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - F64_EPSILON) {
        reverse(solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    } else {
        return false;
    }
}

static bool search_reinsertion(Solution &solut, std::size_t opt,
                               const Data &info) {
    double cost_best = F64_MAX;
    double cost_new;

    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;

    std::size_t I = 0;
    std::size_t J = 0;
    std::size_t POS = 0;

    for (std::size_t i = 1; i < info.dimen - opt + 1; i++) {
        const std::size_t j = opt + i - 1;
        const std::size_t i_prev = i - 1;
        const std::size_t j_next = j + 1;

        const double seq_i_j_w = solut.seq(i, j).w;
        const double seq_i_j_c = solut.seq(i, j).c;
        const double seq_jnext_dimen_w = solut.seq(j_next, info.dimen).w;
        const double seq_jnext_dimen_c = solut.seq(j_next, info.dimen).c;

        for (std::size_t k = 0; k < i_prev; k++) {
            const std::size_t k_next = k + 1;

            cost_concat_1 = solut.seq(0, k).t + info.c(solut[k], solut[i]);
            cost_concat_2 = cost_concat_1 + solut.seq(i, j).t +
                            info.c(solut[j], solut[k_next]);
            cost_concat_3 = cost_concat_2 + solut.seq(k_next, i_prev).t +
                            info.c(solut[i_prev], solut[j_next]);

            cost_new = solut.seq(0, k).c + seq_i_j_w * cost_concat_1 +
                       seq_i_j_c +
                       solut.seq(k_next, i_prev).w * cost_concat_2 +
                       solut.seq(k_next, i_prev).c +
                       seq_jnext_dimen_w * cost_concat_3 + seq_jnext_dimen_c;

            if (cost_new < cost_best) {
                cost_best = cost_new - F64_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (std::size_t k = i + opt; k < info.dimen; k++) {
            const std::size_t k_next = k + 1;

            cost_concat_1 =
                solut.seq(0, i_prev).t + info.c(solut[i_prev], solut[j_next]);
            cost_concat_2 = cost_concat_1 + solut.seq(j_next, k).t +
                            info.c(solut[k], solut[i]);
            cost_concat_3 = cost_concat_2 + solut.seq(i, j).t +
                            info.c(solut[j], solut[k_next]);

            cost_new = solut.seq(0, i_prev).c +
                       solut.seq(j_next, k).w * cost_concat_1 +
                       solut.seq(j_next, k).c +
                       solut.seq(i, j).w * cost_concat_2 + solut.seq(i, j).c +
                       solut.seq(k_next, info.dimen).w * cost_concat_3 +
                       solut.seq(k_next, info.dimen).c;

            if (cost_new < cost_best) {
                cost_best = cost_new - F64_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - F64_EPSILON) {
        reinsert(solut.s, I, J, POS + 1);
        update_subseq_info_matrix(solut, info);
        return true;
    } else {
        return false;
    }
}

static void rvnd(Solution &solut, Data &info) {
    std::vector<Operator> neighbd_list = {Operator::Swap, Operator::TwoOpt,
                                          Operator::Reinsertion,
                                          Operator::OrOpt2, Operator::OrOpt3};

    while (neighbd_list.empty() == false) {
        const std::size_t index = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        const Operator neighbd = neighbd_list[index];

        bool improv_flag;
        switch (neighbd) {
        case Operator::Swap:
            improv_flag = search_swap(solut, info);
            break;
        case Operator::TwoOpt:
            improv_flag = search_two_opt(solut, info);
            break;
        case Operator::Reinsertion:
            improv_flag = search_reinsertion(
                solut, static_cast<std::size_t>(Operator::Reinsertion), info);
            break;
        case Operator::OrOpt2:
            improv_flag = search_reinsertion(
                solut, static_cast<std::size_t>(Operator::OrOpt2), info);
            break;
        case Operator::OrOpt3:
            improv_flag = search_reinsertion(
                solut, static_cast<std::size_t>(Operator::OrOpt3), info);
            break;
        }

        if (improv_flag) {
            neighbd_list = {Operator::Swap, Operator::TwoOpt,
                            Operator::Reinsertion, Operator::OrOpt2,
                            Operator::OrOpt3};
        } else {
            neighbd_list.erase(neighbd_list.begin() + index);
        }
    }
}

static std::vector<std::size_t> perturb(const std::vector<std::size_t> &sl,
                                        Data &info) {
    std::vector<std::size_t> s = sl;
    std::size_t a_start = 1;
    std::size_t a_end = 1;
    std::size_t b_start = 1;
    std::size_t b_end = 1;

    while ((a_start <= b_start && b_start <= a_end) ||
           (b_start <= a_start && a_start <= b_end)) {
        a_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        a_end = a_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

        b_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        b_end = b_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;
    }

    if (a_start < b_start) {
        reinsert(s, b_start, b_end - 1, a_end);
        reinsert(s, a_start, a_end - 1, b_end);
    } else {
        reinsert(s, a_start, a_end - 1, b_end);
        reinsert(s, b_start, b_end - 1, a_end);
    }

    return s;
}

static void gils_rvnd(std::size_t imax, std::size_t iils, const double r[26],
                      Data &info, bool verbose) {
    Solution solut_best{SubseqMatrix(info.dimen),
                        std::vector<std::size_t>(info.dimen, 0), F64_MAX};

    Solution solut_partial{SubseqMatrix(info.dimen),
                           std::vector<std::size_t>(info.dimen, 0), 0.0};

    Solution solut_crnt{SubseqMatrix(info.dimen),
                        std::vector<std::size_t>(info.dimen, 0), 0.0};

    for (std::size_t _i = 0; _i < imax; _i++) {
        const std::size_t r_value = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        const double alpha = r[r_value];

        if (verbose) {
            std::printf("[+] Local Search %zu\n", _i);
        }
        solut_crnt.s = construction(alpha, info);
        if (verbose) {
            std::printf("%s\n", fmt_vec(solut_crnt.s).c_str());
        }

        update_subseq_info_matrix(solut_crnt, info);
        if (verbose) {
            std::printf("\t[+] Constructing Inital Solution.. %s\n",
                        fmt_f64(solut_crnt.cost).c_str());
        }
        solut_partial.s = solut_crnt.s;
        solut_partial.cost = solut_crnt.cost;

        if (verbose) {
            std::printf("\t[+] Looking for the best Neighbor..\n");
        }
        std::size_t iter_ils = 0;
        while (iter_ils < iils) {
            rvnd(solut_crnt, info);
            if (solut_crnt.cost < solut_partial.cost) {
                solut_partial.s = solut_crnt.s;
                solut_partial.cost = solut_crnt.cost;
                iter_ils = 0;
            }

            solut_crnt.s = perturb(solut_partial.s, info);
            update_subseq_info_matrix(solut_crnt, info);
            iter_ils += 1;
        }

        if (solut_partial.cost < solut_best.cost) {
            solut_best.s = solut_partial.s;
            solut_best.cost = solut_partial.cost;
        }

        if (verbose) {
            std::printf("\tCurrent Best Cost %s\n",
                        fmt_f64(solut_best.cost).c_str());
        }
    }

    std::printf("%s\n", fmt_vec(solut_best.s).c_str());
    std::printf("COST: %s\n", fmt_f64(solut_best.cost).c_str());
}

int main(int argc, char **argv) {
    const std::string filename = "../distance_matrix";

    bool verbose = false;
    for (int i = 0; i < argc; i++) {
        const std::string arg(argv[i]);
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
    }

    data::LoadResult loaded = data::load(filename);

    Data info{std::move(loaded.c), loaded.dimension, std::move(loaded.rnd), 0};

    const std::size_t imax = 10;
    const std::size_t iils = (loaded.dimension < 100) ? loaded.dimension : 100;

    const double r[26] = {
        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11,
        0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23,
        0.24, 0.25,
    };

    timespec start;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    gils_rvnd(imax, iils, r, info, verbose);

    timespec end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    const double cpu_time = static_cast<double>(end.tv_sec - start.tv_sec) +
                            static_cast<double>(end.tv_nsec - start.tv_nsec) /
                                1e9;
    std::printf("TIME: %s\n", fmt_f64_debug(cpu_time).c_str());

    return 0;
}
