#include "Data.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>

int readData(std::vector<std::vector<double>> &c, std::vector<int> &rnd, int &rnd_size) {
    const char *filename = "../distance_matrix";
    FILE *fh = fopen(filename, "rb");
    if (fh == NULL) {
        fprintf(stderr, "%s: No such file or directory\n", filename);
        exit(1);
    }

    std::string contents;
    {
        char buf[65536];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), fh)) > 0) {
            contents.append(buf, n);
        }
    }
    fclose(fh);

    int index = 1;
    int dimension = 0;
    // 1-based list of the non-empty lines, like string.gmatch(contents, "[^\n]+")
    std::vector<std::string> lines;
    lines.push_back(std::string());
    {
        size_t pos = 0;
        while (pos < contents.size()) {
            size_t nl = contents.find('\n', pos);
            if (nl == std::string::npos) nl = contents.size();
            if (nl > pos) {
                lines.push_back(contents.substr(pos, nl - pos));
            }
            pos = nl + 1;
        }
    }

    dimension = atoi(lines[index].c_str());
    index = index + 1;

    c.assign(dimension + 1, std::vector<double>(dimension + 1, 0.0));

    for (int i = 1; i <= dimension - 1; i++) {
        int j = i + 1;
        std::istringstream iss(lines[index]);
        std::string cost;
        while (iss >> cost) {
            c[i][j] = atof(cost.c_str());
            c[j][i] = atof(cost.c_str());
            j = j + 1;
        }
        index = index + 1;
    }

    index = index + 3;

    rnd_size = atoi(lines[index].c_str());
    printf("%s\t%d\n", lines[index].c_str(), (int)(rnd.size() ? rnd.size() - 1 : 0));

    rnd.assign(1, 0);
    for (int i = 1; i <= rnd_size; i++) {
        rnd.push_back(atoi(lines[index + i].c_str()));
    }
    printf("%s\n", lines[index + rnd_size].c_str());

    return dimension;
}
