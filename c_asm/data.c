#include "data.h"

static
Real ** matrix_allocate(size_t sz) {
    Real ** ptr = (Real **) calloc(sz, sizeof(Real*));

    if (ptr == NULL) printf("qebro\n");
    for ( int i = 0; i < sz; i++ ) {
        //ptr[i] = new double [sz];
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
        //string line;
        //char buff[3000];

        fscanf(file,"%d", &dimension);
        printf("dimension:  %d\n", dimension);
        //getline(file, line);
        //dimension = stoi(line);
        *matrix = matrix_allocate((size_t) dimension);
        //cout << "dimension " << dimension << endl;

        _Bool flag = false;
        char * read_fmt;
        if      (sizeof(Real) == sizeof(float))  read_fmt = "%f";
        else if (sizeof(Real) == sizeof(double)) read_fmt = "%lf";
        for (int i = 0; i < dimension-1; i++) {
            //cout << line << "\n";
            //int j = i + 1;
            for (int j = i+1; j < dimension; j++) {
                Real value;
                fscanf(file, read_fmt, &value);
                    
                (*matrix)[i][j] = value;
                (*matrix)[j][i] = value;
                //printf("%.2lf ", value);
                //cout << line << endl;
            }
            //printf("\n");

        }


        char buff[100];
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        fgets(buff, 100, file);
        //printf("%s", buff);
      //fprintf(file, NULL, NULL);
      //fprintf(file, NULL, NULL);
      //fprintf(file, NULL, NULL);


        int rnd_size;
        fscanf(file, "%d", &rnd_size);
        //printf("%d\n", rnd_size);
        
        *rnd = (int *) calloc(rnd_size, sizeof(int));

        //if (rnd == NULL) printf("eita\n");
        for (int i = 0; i < rnd_size; i++) {
            int value;
            fscanf(file, "%d", &value);
            (*rnd)[i] = value;
            //printf("%d\n", (*rnd)[i]);
        }
        //exit(0);

    }

    fclose(file);
    return dimension;
}
