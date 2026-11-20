#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <sys/resource.h>

static double cputime() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    double u = double(r.ru_utime.tv_sec) + double(r.ru_utime.tv_usec) / 1e6;
    double s = double(r.ru_stime.tv_sec) + double(r.ru_stime.tv_usec) / 1e6;
    return u + s;
}

/* formats a float the same way Go's fmt.Println does (strconv 'g', prec -1) */
static std::string fmt_float(double v) {
    if (std::isnan(v)) {
        return "NaN";
    }
    if (std::isinf(v)) {
        return v > 0 ? "+Inf" : "-Inf";
    }
    if (v == 0.0) {
        return std::signbit(v) ? "-0" : "0";
    }

    char buf[64];
    int p = 0;
    for (; p < 17; p++) {
        snprintf(buf, sizeof(buf), "%.*e", p, v);
        if (strtod(buf, nullptr) == v) {
            break;
        }
    }

    /* buf holds [-]d[.ddd]e[+-]dd : split it into digits and exponent */
    std::string str(buf);
    bool neg = false;
    size_t pos = 0;
    if (str[0] == '-') {
        neg = true;
        pos = 1;
    }

    std::string digits;
    size_t e = str.find('e');
    for (size_t i = pos; i < e; i++) {
        if (str[i] != '.') {
            digits.push_back(str[i]);
        }
    }
    while (digits.size() > 1 && digits.back() == '0') {
        digits.pop_back();
    }
    int exp = atoi(str.c_str() + e + 1);

    std::string out;
    if (exp < -4 || exp >= 6) {
        snprintf(buf, sizeof(buf), "%.*e", int(digits.size()) - 1, v);
        return std::string(buf);
    }

    int nd = int(digits.size());
    int dp = exp + 1;
    if (dp <= 0) {
        out = "0.";
        out.append(size_t(-dp), '0');
        out += digits;
    } else if (dp >= nd) {
        out = digits;
        out.append(size_t(dp - nd), '0');
    } else {
        out = digits.substr(0, size_t(dp)) + "." + digits.substr(size_t(dp));
    }

    return neg ? "-" + out : out;
}

/* formats an int slice the same way Go's fmt.Println does */
static std::string fmt_slice(const std::vector<int> &v) {
    std::string out = "[";
    for (size_t i = 0; i < v.size(); i++) {
        if (i > 0) {
            out += " ";
        }
        out += std::to_string(v[i]);
    }
    out += "]";
    return out;
}

/* Mimics Go's fmt.Fscanf verbs on a file: a newline read while scanning ends
 * the current scan (the newline is consumed), so a scan never crosses a line. */
struct tScanner {
    FILE *f;

    /* "%d" ; leaves *out untouched and returns false when the scan fails */
    bool scan_int(int *out) {
        int ch;
        for (;;) {
            ch = getc(f);
            if (ch == EOF) {
                return false;
            }
            if (ch == '\n') {
                return false;
            }
            if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\v' || ch == '\f') {
                continue;
            }
            break;
        }

        bool neg = false;
        if (ch == '-' || ch == '+') {
            neg = (ch == '-');
            ch = getc(f);
        }

        long v = 0;
        bool any = false;
        while (ch >= '0' && ch <= '9') {
            any = true;
            v = v * 10 + (ch - '0');
            ch = getc(f);
        }

        if (!any) {
            return false;
        }
        if (ch != EOF && ch != '\n') {
            ungetc(ch, f);
        }

        *out = int(neg ? -v : v);
        return true;
    }

    /* "%s\n" ; leaves *out untouched and returns false when the scan fails */
    bool scan_token(std::string *out) {
        int ch;
        for (;;) {
            ch = getc(f);
            if (ch == EOF) {
                return false;
            }
            if (ch == '\n') {
                return false;
            }
            if (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\v' || ch == '\f') {
                continue;
            }
            break;
        }

        std::string t;
        while (ch != EOF && ch != '\n' && ch != ' ' && ch != '\t' && ch != '\r' && ch != '\v' && ch != '\f') {
            t.push_back(char(ch));
            ch = getc(f);
        }
        if (ch != EOF && ch != '\n') {
            ungetc(ch, f);
        }

        *out = t;
        return true;
    }
};

struct tSubseq {
    double T;
    double W;
    double C;
};

struct tData {
    std::vector<std::vector<double>> c;
    int dimen;
    std::vector<int> rnd;
    int rnd_index;
};

struct tSolution {
    std::vector<int> s;
    std::vector<tSubseq> seq;
    int dim;
    double cost;

    tSubseq *at(int i, int j) { return &seq[i * dim + j]; }
};

enum {
    SWAP = 0,
    REINSERTION = 1,
    OR_OPT_2 = 2,
    OR_OPT_3 = 3,
    TWO_OPT = 4,
};

static void read_data(bool verbose, int &dimen_out, std::vector<std::vector<double>> &cost_out,
                      std::vector<int> &rnd_out) {
    std::vector<std::vector<double>> cost;
    int dimen = 0;
    int c = 0;
    int rnd_size = 0;

    FILE *fp = fopen("../distance_matrix", "r");
    if (fp == nullptr) {
        perror("../distance_matrix");
        exit(1);
    }
    tScanner file{fp};

    file.scan_int(&dimen);

    if (verbose) {
        std::cout << dimen << std::endl;
    }

    cost.resize(dimen);

    for (int i = 0; i < dimen; i++) {
        cost[i].resize(dimen);
        for (int j = i + 1; j < dimen; j++) {
            file.scan_int(&c);

            cost[i][j] = double(c);
            if (verbose) {
                printf("%f ", cost[i][j]);
            }
        }
        // o ultimo valor da linha é lido duas vezes para que a função inicie a leitura da próxima linha do arquivo.
        file.scan_int(&c);
        if (verbose) {
            std::cout << std::endl;
        }
    }

    for (int i = 0; i < dimen; i++) {
        for (int j = i + 1; j < dimen; j++) {
            cost[j][i] = cost[i][j];
        }
    }

    std::string buff;
    file.scan_token(&buff);
    file.scan_token(&buff);
    file.scan_token(&buff);
    if (verbose) {
        std::cout << buff << std::endl;
    }
    file.scan_int(&rnd_size);
    if (verbose) {
        std::cout << rnd_size << std::endl;
    }

    std::vector<int> rnd(rnd_size);

    for (int i = 0; i < rnd_size; i++) {
        file.scan_int(&rnd[i]);
    }

    fclose(fp);

    dimen_out = dimen;
    cost_out = cost;
    rnd_out = rnd;
}

static bool feasible_(const std::vector<int> &s, const tData &data) {
    std::vector<bool> is(data.dimen);

    for (int i = 0; i < data.dimen; i++) {
        is[i] = false;
    }

    for (int i = 0; i < data.dimen; i++) {
        is[s[i]] = true;
    }

    for (int i = 0; i < data.dimen; i++) {
        if (is[i] == false) {
            return false;
        }
    }

    return true;
}

static void remove_at(std::vector<int> &arr, int i) { arr.erase(arr.begin() + i); }

static int partition(std::vector<int> &arr, int left, int right, tData *data, int r) {
    int pivot = arr[right];
    int i = left - 1;
    for (int j = left; j < right; j++) {
        if (data->c[r][arr[j]] < data->c[r][pivot]) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[right]);
    return i + 1;
}

static void quicksort(std::vector<int> &arr, int left, int right, tData *data, int r) {
    if (left < right) {
        int pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

static void sort_(std::vector<int> &arr, int r, tData *data) {
    quicksort(arr, 0, int(arr.size()) - 1, data, r);
}

static std::vector<int> construction(double alpha, tData *data) {
    (void)alpha;

    std::vector<int> s(1);
    s[0] = 0;

    std::vector<int> cL(data->dimen - 1);
    for (int i = 0; i < data->dimen - 1; i++) {
        cL[i] = i + 1;
    }

    int r = 0;
    while (cL.size() > 0) {
        sort_(cL, r, data);

        int index = data->rnd[data->rnd_index];
        data->rnd_index++;

        int c = cL[index];
        r = c;

        remove_at(cL, index);
        s.push_back(c);
    }

    s.push_back(0);

    return s;
}

static tSolution NewSolution(const tData &data) {
    tSolution solut;

    solut.s.assign(data.dimen + 1, 0);
    solut.dim = data.dimen + 1;
    solut.seq.assign(size_t(solut.dim) * size_t(solut.dim), tSubseq{0.0, 0.0, 0.0});
    solut.cost = 0.0;

    return solut;
}

static void update_subseq_info_matrix(tSolution *solut, tData *data) {

    for (int i = 0; i < data->dimen + 1; i++) {
        int k = 1 - i;

        if (i == 0) {
            k--;
        }

        solut->at(i, i)->T = 0.0;
        solut->at(i, i)->C = 0.0;
        if (i == 0) {
            solut->at(i, i)->W = 0.0;
        } else {
            solut->at(i, i)->W = 1.0;
        }

        for (int j = i + 1; j < data->dimen + 1; j++) {
            int j_prev = j - 1;

            solut->at(i, j)->T = data->c[solut->s[j_prev]][solut->s[j]] + solut->at(i, j_prev)->T;

            solut->at(i, j)->C = solut->at(i, j)->T + solut->at(i, j_prev)->C;

            solut->at(i, j)->W = double(j + k);
        }
    }

    solut->cost = solut->at(0, data->dimen)->C;
}

static void swap_(tSolution *solut, int i, int j) {
    int tmp = solut->s[i];
    solut->s[i] = solut->s[j];
    solut->s[j] = tmp;
}

static bool search_swap(tSolution *solut, tData *data) {

    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;
    double cost_concat_4;

    double cost_best = std::numeric_limits<double>::max();
    double cost_new;
    int I = 0;
    int J = 0;

    for (int i = 1; i < data->dimen - 1; i++) {
        int i_prev = i - 1;
        int i_next = i + 1;

        cost_concat_1 = solut->at(0, i_prev)->T + data->c[solut->s[i_prev]][solut->s[i_next]];
        cost_concat_2 = cost_concat_1 + solut->at(i, i_next)->T + data->c[solut->s[i]][solut->s[i_next + 1]];

        cost_new = solut->at(0, i_prev)->C +
                   solut->at(i, i_next)->W * cost_concat_1 + data->c[solut->s[i_next]][solut->s[i]] +
                   solut->at(i_next + 1, data->dimen)->W * cost_concat_2 + solut->at(i_next + 1, data->dimen)->C;

        if (cost_new < cost_best) {
            cost_best = cost_new;
            I = i;
            J = i_next;
        }

        for (int j = i_next + 1; j < data->dimen; j++) {
            int j_next = j + 1;
            int j_prev = j - 1;

            cost_concat_1 = solut->at(0, i_prev)->T + data->c[solut->s[i_prev]][solut->s[j]];
            cost_concat_2 = cost_concat_1 + data->c[solut->s[j]][solut->s[i_next]];
            cost_concat_3 = cost_concat_2 + solut->at(i_next, j_prev)->T + data->c[solut->s[j_prev]][solut->s[i]];
            cost_concat_4 = cost_concat_3 + data->c[solut->s[i]][solut->s[j_next]];

            cost_new = solut->at(0, i_prev)->C +
                       cost_concat_1 +
                       solut->at(i_next, j_prev)->W * cost_concat_2 + solut->at(i_next, j_prev)->C +
                       cost_concat_3 +
                       solut->at(j_next, data->dimen)->W * cost_concat_4 + solut->at(j_next, data->dimen)->C;

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->at(0, data->dimen)->C) {
        swap_(solut, I, J);

        update_subseq_info_matrix(solut, data);
        return true;
    }

    return false;
}

static void reverse(tSolution *solut, int i, int j) {
    int f = i;
    int l = j;

    int m = int((i + j) / 2);

    while (f <= m) {
        swap_(solut, f, l);
        f++;
        l--;
    }
}

static bool search_two_opt(tSolution *solut, tData *data) {
    double cost_new;
    double cost_best = std::numeric_limits<double>::max();

    double cost_concat_1;
    double cost_concat_2;

    int I = 0;
    int J = 0;

    for (int i = 1; i < data->dimen - 1; i++) {
        int i_prev = i - 1;
        double rev_seq_cost = solut->at(i, i + 1)->T;

        for (int j = i + 2; j < data->dimen; j++) {
            int j_next = j + 1;

            rev_seq_cost += data->c[solut->s[j - 1]][solut->s[j]] * (solut->at(i, j)->W - 1.0);

            cost_concat_1 = solut->at(0, i_prev)->T + data->c[solut->s[j]][solut->s[i_prev]];
            cost_concat_2 = cost_concat_1 + solut->at(i, j)->T + data->c[solut->s[j_next]][solut->s[i]];

            cost_new = solut->at(0, i_prev)->C +
                       solut->at(i, j)->W * cost_concat_1 + rev_seq_cost +
                       solut->at(j_next, data->dimen)->W * cost_concat_2 + solut->at(j_next, data->dimen)->C;

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut->cost) {
        reverse(solut, I, J);

        update_subseq_info_matrix(solut, data);

        return true;
    }

    return false;
}

static void reinsert(tSolution *solut, int i, int j, int pos) {
    int sz = j - i + 1;
    std::vector<int> sub(sz);

    memcpy(sub.data(), &solut->s[i], size_t(sz) * sizeof(int));

    if (pos < i) {
        memmove(&solut->s[pos + sz], &solut->s[pos], size_t(i - pos) * sizeof(int));
        memcpy(&solut->s[pos], sub.data(), size_t(sz) * sizeof(int));

    } else {
        memmove(&solut->s[i], &solut->s[j + 1], size_t(pos - j - 1) * sizeof(int));
        memcpy(&solut->s[pos - sz], sub.data(), size_t(sz) * sizeof(int));
    }
}

static bool search_reinsertion(tSolution *solut, tData *data, int opt) {
    double cost_best = std::numeric_limits<double>::max();
    double cost_new;

    double cost_concat_1;
    double cost_concat_2;
    double cost_concat_3;

    int I = 0;
    int J = 0;
    int POS = 0;

    for (int i = 1; i < data->dimen - opt + 1; i++) {
        int j = opt + i - 1;
        int i_prev = i - 1;
        int j_next = j + 1;

        for (int k = 0; k < i_prev; k++) {
            int k_next = k + 1;

            cost_concat_1 = solut->at(0, k)->T + data->c[solut->s[k]][solut->s[i]];
            cost_concat_2 = cost_concat_1 + solut->at(i, j)->T + data->c[solut->s[j]][solut->s[k_next]];
            cost_concat_3 = cost_concat_2 + solut->at(k_next, i_prev)->T + data->c[solut->s[i_prev]][solut->s[j_next]];

            cost_new = solut->at(0, k)->C +                                                    /*        1st subseq */
                       solut->at(i, j)->W * cost_concat_1 + solut->at(i, j)->C +               /* concat 2nd subseq (reinserted seq) */
                       solut->at(k_next, i_prev)->W * cost_concat_2 + solut->at(k_next, i_prev)->C + /* concat 3rd subseq */
                       solut->at(j_next, data->dimen)->W * cost_concat_3 + solut->at(j_next, data->dimen)->C; /* concat 4th subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (int k = i + opt; k < data->dimen; k++) {
            int k_next = k + 1;

            cost_concat_1 = solut->at(0, i_prev)->T + data->c[solut->s[i_prev]][solut->s[j_next]];
            cost_concat_2 = cost_concat_1 + solut->at(j_next, k)->T + data->c[solut->s[k]][solut->s[i]];
            cost_concat_3 = cost_concat_2 + solut->at(i, j)->T + data->c[solut->s[j]][solut->s[k_next]];

            cost_new = solut->at(0, i_prev)->C +                                               /*      1st subseq */
                       solut->at(j_next, k)->W * cost_concat_1 + solut->at(j_next, k)->C +     /* concat 2nd subseq */
                       solut->at(i, j)->W * cost_concat_2 + solut->at(i, j)->C +               /* concat 3rd subseq (reinserted seq) */
                       solut->at(k_next, data->dimen)->W * cost_concat_3 + solut->at(k_next, data->dimen)->C; /* concat 4th subseq */

            if (cost_new < cost_best) {
                cost_best = cost_new;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut->cost) {
        reinsert(solut, I, J, POS + 1);

        update_subseq_info_matrix(solut, data);

        return true;
    }

    return false;
}

static void RVND(tSolution *solut, tData *data) {
    const std::vector<int> n_list_b = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3};
    std::vector<int> n_list = n_list_b;

    while (n_list.size() > 0) {
        int index = data->rnd[data->rnd_index];
        data->rnd_index++;

        bool improve = false;
        switch (n_list[index]) {
        case REINSERTION:
            improve = search_reinsertion(solut, data, REINSERTION);
            break;
        case OR_OPT_2:
            improve = search_reinsertion(solut, data, OR_OPT_2);
            break;
        case OR_OPT_3:
            improve = search_reinsertion(solut, data, OR_OPT_3);
            break;
        case TWO_OPT:
            improve = search_two_opt(solut, data);
            break;
        case SWAP:
            improve = search_swap(solut, data);
            break;
        }

        if (improve == true) {
            n_list = n_list_b;

        } else {
            remove_at(n_list, index);
        }
    }
}

static std::vector<int> perturb(const std::vector<int> &sl, tData *data) {
    std::vector<int> s(data->dimen + 1);
    std::copy(sl.begin(), sl.begin() + std::min(sl.size(), s.size()), s.begin());

    int A_start = 1;
    int A_end = 1;
    int B_start = 1;
    int B_end = 1;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data->rnd[data->rnd_index];
        data->rnd_index++;
        A_end = A_start + data->rnd[data->rnd_index];
        data->rnd_index++;

        B_start = data->rnd[data->rnd_index];
        data->rnd_index++;
        B_end = B_start + data->rnd[data->rnd_index];
        data->rnd_index++;
    }

    auto reinsert_s = [&](std::vector<int> *s_, int i, int j, int pos) {
        int sz = j - i + 1;
        std::vector<int> sub(j - i + 1);

        memcpy(sub.data(), &(*s_)[i], size_t(sz) * sizeof(int));

        if (pos < i) {
            memmove(&(*s_)[pos + sz], &(*s_)[pos], size_t(i - pos) * sizeof(int));
            memcpy(&(*s_)[pos], sub.data(), size_t(sz) * sizeof(int));

        } else {
            memmove(&(*s_)[i], &(*s_)[j + 1], size_t(pos - j - 1) * sizeof(int));
            memcpy(&(*s_)[pos - sz], sub.data(), size_t(sz) * sizeof(int));
        }

        if (feasible_(*s_, *data) == false) {
            exit(0);
        }
    };

    if (A_start < B_start) {
        reinsert_s(&s, B_start, B_end - 1, A_end);
        reinsert_s(&s, A_start, A_end - 1, B_end);

    } else {
        reinsert_s(&s, A_start, A_end - 1, B_end);
        reinsert_s(&s, B_start, B_end - 1, A_end);
    }

    return s;
}

static void GILS_RVND(int Imax, int Iils, const std::array<double, 26> &R, tData data, bool verbose) {

    tSolution solut_crnt = NewSolution(data);
    tSolution solut_partial = NewSolution(data);
    tSolution solut_best = NewSolution(data);
    solut_best.cost = std::numeric_limits<double>::infinity();

    for (int i = 0; i < Imax; i++) {

        if (verbose) {
            printf("[+] Local Search %d\n", i);
        }

        int index = data.rnd[data.rnd_index];
        data.rnd_index++;

        solut_crnt.s = construction(R[index], &data);
        update_subseq_info_matrix(&solut_crnt, &data);
        if (verbose) {
            printf("\t[+] Constructing Inital Solution.. %.2f\n", solut_crnt.cost);
            std::cout << "\t " << fmt_slice(solut_crnt.s) << std::endl;
        }

        solut_partial.s = solut_crnt.s;
        solut_partial.cost = solut_crnt.cost;

        if (verbose) {
            std::cout << "\t[+] Looking for the best Neighbor.." << std::endl;
        }
        int iterILS = 0;
        while (iterILS < Iils) {
            RVND(&solut_crnt, &data);
            if (solut_crnt.cost < solut_partial.cost) {
                solut_partial.s = solut_crnt.s;
                solut_partial.cost = solut_crnt.cost;
                iterILS = 0;
            }

            solut_crnt.s = perturb(solut_partial.s, &data);
            update_subseq_info_matrix(&solut_crnt, &data);
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost) {
            solut_best.cost = solut_partial.cost;
            solut_best.s = solut_partial.s;
        }
    }

    std::cout << "COST: " << " " << fmt_float(solut_best.cost) << std::endl;
    std::cout << "SOLUTION: " << " " << fmt_slice(solut_best.s) << std::endl;
}

int main(int argc, char **argv) {
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
    }

    if (verbose) {
        std::cout << "Hello Vourld!" << std::endl;
    }

    tData data;
    data.rnd_index = 0;

    read_data(verbose, data.dimen, data.c, data.rnd);

    std::array<double, 26> R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                                0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};

    int Imax = 10;
    int Iils = 100;
    if (data.dimen < 100) {
        Iils = data.dimen;
    }

    double cpuStart = cputime();
    auto start = std::chrono::steady_clock::now();

    GILS_RVND(Imax, Iils, R, data, verbose);

    double cpuElapsed = cputime() - cpuStart;
    auto t = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t - start).count();
    std::cout << "TIME: " << fmt_float(cpuElapsed) << std::endl;
    std::cout << "wall clock (s): " << fmt_float(elapsed) << std::endl;

    return 0;
}
