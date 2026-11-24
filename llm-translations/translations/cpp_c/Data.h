#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static
double ** matrix_allocate(size_t sz) {
    double ** ptr = (double **) malloc(sz * sizeof(double *));

    for ( size_t i = 0; i < sz; i++ ) {
        ptr[i] = (double *) calloc(sz, sizeof(double));
    }
    return ptr;
}

static
int loadData(double *** matrix, int ** rnd, int * rnd_len) {
    int dimension = 0;
    FILE * file;

    *rnd = NULL;
    *rnd_len = 0;

    file = fopen("../distance_matrix", "r");

    if (file != NULL) {
        char * line = NULL;
        size_t cap = 0;
        ssize_t len;

        if ((len = getline(&line, &cap, file)) != -1) {
            while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
                line[--len] = '\0';

            dimension = atoi(line);
            *matrix = matrix_allocate((size_t) dimension);

            int flag = 0;
            int rnd_size = -1;
            int rnd_count = 0;
            int i = 0;
            while ((len = getline(&line, &cap, file)) != -1) {
                while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
                    line[--len] = '\0';

                int j = i + 1;
                char * cur = line;
                char * pos;
                while ((pos = strchr(cur, ' ')) != NULL) {
                    *pos = '\0';
                    double value = atof(cur);
                    cur = pos + 1;
                    (*matrix)[i][j] = value;
                    (*matrix)[j][i] = value;
                    j++;
                }

                i++;

                if (flag) {

                    if (rnd_size < 0) {
                        rnd_size = atoi(cur);
                        *rnd = (int *) malloc((size_t) rnd_size * sizeof(int));
                    } else {
                        (*rnd)[rnd_count++] = atoi(cur);
                    }
                }

                if (strcmp(cur, "RND") == 0) flag = 1;
            }

            *rnd_len = rnd_count;
        }

        free(line);
        fclose(file);
    }

    return dimension;
}

#endif
