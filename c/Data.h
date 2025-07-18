
#define TRUE 1
#define FALSE 0

double ** matrix_allocate(size_t sz) {
    double ** ptr = (double **) calloc(sz, sizeof(double*));

    if (ptr == NULL) printf("qebro\n");
    for ( int i = 0; i < sz; i++ ) {
        ptr[i] = (double *) calloc(sz, sizeof(double));
        if (ptr[i] == NULL) printf("qebro\n");
    }

    return ptr;
}

int loadData(double *** matrix, int ** rnd) {
    int dimension;
    FILE * file;

    file = fopen("../distance_matrix", "r");

    if (file != NULL) {
        fscanf(file,"%d", &dimension);
        printf("dimension:  %d\n", dimension);
        *matrix = matrix_allocate((size_t) dimension);

        char flag = FALSE;
        for (int i = 0; i < dimension-1; i++) {
            for (int j = i+1; j < dimension; j++) {
                double value;
                fscanf(file, "%lf", &value);
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
