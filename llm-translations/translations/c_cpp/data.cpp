#include "data.hpp"

static
std::vector<std::vector<Real>> matrix_allocate(size_t sz) {
    return std::vector<std::vector<Real>>(sz, std::vector<Real>(sz, 0.0));
}

int loadData(std::vector<std::vector<Real>> & matrix, std::vector<int> & rnd) {
    int dimension = 0;
    FILE * file;

    file = fopen("../distance_matrix", "r");

    if (file != NULL) {
        fscanf(file, "%d", &dimension);
        printf("dimension:  %d\n", dimension);
        matrix = matrix_allocate((size_t) dimension);

        const char * read_fmt = "%lf";
        if      (sizeof(Real) == sizeof(float))  read_fmt = "%f";
        else if (sizeof(Real) == sizeof(double)) read_fmt = "%lf";
        for (int i = 0; i < dimension-1; i++) {
            for (int j = i+1; j < dimension; j++) {
                Real value;
                fscanf(file, read_fmt, &value);

                matrix[i][j] = value;
                matrix[j][i] = value;
            }
        }


        char buff[100];
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);


        int rnd_size;
        fscanf(file, "%d", &rnd_size);

        rnd.assign(rnd_size, 0);

        for (int i = 0; i < rnd_size; i++) {
            int value;
            fscanf(file, "%d", &value);
            rnd[i] = value;
        }

    }

    fclose(file);
    return dimension;
}
