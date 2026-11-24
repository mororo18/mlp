#include "data.h"

static
Real ** matrix_allocate(size_t sz) {
    Real ** ptr = (Real **) calloc(sz, sizeof(Real*));

    if (ptr == NULL) printf("qebro\n");
    for ( int i = 0; i < sz; i++ ) {
        ptr[i] = (Real *) calloc(sz, sizeof(Real));
        if (ptr[i] == NULL) printf("qebro\n");
    }

    return ptr;
}

int loadData(Real *** matrix, int ** rnd) {
    int dimension;
    FILE * file;

    file = fopen("../distance_matrix", "r");

    if (file != NULL) {
        fscanf(file,"%d", &dimension);
        printf("dimension:  %d\n", dimension);
        *matrix = matrix_allocate((size_t) dimension);

        _Bool flag = false;
        char * read_fmt;
        if      (sizeof(Real) == sizeof(float))  read_fmt = "%f";
        else if (sizeof(Real) == sizeof(double)) read_fmt = "%lf";
        for (int i = 0; i < dimension-1; i++) {
            for (int j = i+1; j < dimension; j++) {
                Real value;
                fscanf(file, read_fmt, &value);
                    
                (*matrix)[i][j] = value;
                (*matrix)[j][i] = value;
            }
        }


        char buff[100];
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);


        int rnd_size;
        fscanf(file, "%d", &rnd_size);
        
        *rnd = (int *) calloc(rnd_size, sizeof(int));

        for (int i = 0; i < rnd_size; i++) {
            int value;
            fscanf(file, "%d", &value);
            (*rnd)[i] = value;
        }

    }

    fclose(file);
    return dimension;
}
