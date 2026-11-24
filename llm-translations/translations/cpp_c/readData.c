#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "readData.h"

double CalcDistEuc ( double *X, double *Y, int I, int J );
double CalcDistAtt ( double *X, double *Y, int I, int J );
void CalcLatLong ( double *X, double *Y, int n, double *latit, double* longit );
double CalcDistGeo ( double *latit, double *longit, int I, int J );

void readData( int argc, char** argv, int* Dimension, double ***Mdist )
{

     if (argc < 2) {
         printf("\nFaltando parametros\n");
         printf(" ./exec [Instancia] \n");
         exit(1);
     }

     if (argc > 2) {
          printf("\nMuitos parametros\n");
          printf(" ./exec [Instancia] \n");
         exit(1);
     }

    int N;
    char arquivo[256], ewt[256];

     char *instancia;
     instancia = argv[1];


    FILE *in = fopen( instancia, "r");

	if (!in) {
		printf("arquivo nao pode ser aberto\n");
		exit(1);
    }

    arquivo[0] = '\0';

    while ( strcmp(arquivo, "DIMENSION:") != 0 && strcmp(arquivo, "DIMENSION" ) != 0 ) {
        if (fscanf(in, "%255s", arquivo) != 1) break;
    }

    if ( strcmp(arquivo, "DIMENSION" ) == 0 )  fscanf(in, "%255s", arquivo);

    fscanf(in, "%d", &N);

    while ( strcmp(arquivo, "EDGE_WEIGHT_TYPE:") != 0 && strcmp(arquivo, "EDGE_WEIGHT_TYPE" ) != 0 ) {
        if (fscanf(in, "%255s", arquivo) != 1) break;
    }
    if ( strcmp(arquivo, "EDGE_WEIGHT_TYPE" ) == 0 )  fscanf(in, "%255s", arquivo);

    fscanf(in, "%255s", ewt);

    double *x = (double *) malloc((size_t)(N+1) * sizeof(double));
    double *y = (double *) malloc((size_t)(N+1) * sizeof(double));

    // Alocar matriz 2D
    double **dist = (double **) malloc((size_t)(N+1) * sizeof(double *));

    for ( int i = 0; i < N+1; i++ ) {
        dist [i] = (double *) malloc((size_t)(N+1) * sizeof(double));
    }

    if ( strcmp(ewt, "EXPLICIT") == 0 ) {

        while ( strcmp(arquivo, "EDGE_WEIGHT_FORMAT:") != 0 && strcmp(arquivo, "EDGE_WEIGHT_FORMAT" ) != 0 ) {
            if (fscanf(in, "%255s", arquivo) != 1) break;
        }

        char ewf[256];
        if ( strcmp(arquivo, "EDGE_WEIGHT_FORMAT" ) == 0 )  fscanf(in, "%255s", arquivo);
        fscanf(in, "%255s", ewf);

        if ( strcmp(ewf, "FUNCTION") == 0 ) {
            printf("FUNCTION - Nao suportado!\n"); }

        else if ( strcmp(ewf, "FULL_MATRIX") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = 1; j < N+1; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                }
            }
        }

        else if ( strcmp(ewf, "UPPER_ROW") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N; i++ ) {
                for ( int j = i+1; j < N+1; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( strcmp(ewf, "LOWER_ROW") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int i = 2; i < N+1; i++ ) {
                for ( int j = 1; j < i; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }
        }

        else if ( strcmp(ewf, "UPPER_DIAG_ROW") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = i; j < N+1; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( strcmp(ewf, "LOWER_DIAG_ROW") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = 1; j <= i; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( strcmp(ewf, "UPPER_COL") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int j = 2; j < N+1; j++ ) {
                for ( int i = 1; i < j; i++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( strcmp(ewf, "LOWER_COL") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N; j++ ) {
                for ( int i = j+1; i < N+1; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( strcmp(ewf, "UPPER_DIAG_COL") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N+1; j++ ) {
                for ( int i = 1; i <= j; i++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( strcmp(ewf, "LOWER_DIAG_COL") == 0 ) {

            while ( strcmp(arquivo, "EDGE_WEIGHT_SECTION") != 0 ) {
                if (fscanf(in, "%255s", arquivo) != 1) break;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N+1; j++ ) {
                for ( int i = j; i < N+1; j++ ) {
                    fscanf(in, "%lf", &dist[i][j]);
                    dist[j][i] = dist[i][j];
                }
            }
        }

    }

    else if ( strcmp(ewt, "EUC_2D") == 0 ) {

        while ( strcmp(arquivo, "NODE_COORD_SECTION") != 0 ) {
            if (fscanf(in, "%255s", arquivo) != 1) break;
        }
        // ler coordenadas
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            fscanf(in, "%d %lf %lf", &tempCity, &x[i], &y[i]);
        }

        // Calcular Matriz Distancia (Euclidiana)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = floor ( CalcDistEuc ( x, y, i, j ) + 0.5 );
            }
        }
    }

    else if ( strcmp(ewt, "EUD_3D") == 0 ) {
        printf("EUC_3D - Nao suportado!\n"); }

    else if ( strcmp(ewt, "MAX_2D") == 0 ) {
        printf("MAX_2D - Nao suportado!\n"); }

    else if ( strcmp(ewt, "MAX_3D") == 0 ) {
        printf("MAX_3D - Nao Suportado!\n"); }

    else if ( strcmp(ewt, "MAN_2D") == 0 ) {
        printf("MAN_2D - Nao suportado!\n"); }

    else if ( strcmp(ewt, "MAN_3D") == 0 ) {
        printf("MAN_3D - Nao Suportado!\n"); }

    else if ( strcmp(ewt, "CEIL_2D") == 0 ) {

        while ( strcmp(arquivo, "NODE_COORD_SECTION") != 0 ) {
            if (fscanf(in, "%255s", arquivo) != 1) break;
        }
        // ler coordenadas
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            fscanf(in, "%d %lf %lf", &tempCity, &x[i], &y[i]);
        }

        // Calcular Matriz Distancia (Euclidiana)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = ceil ( CalcDistEuc ( x, y, i, j ) );
            }
        }
    }

    else if ( strcmp(ewt, "GEO") == 0 ) {

        while ( strcmp(arquivo, "NODE_COORD_SECTION") != 0 ) {
            if (fscanf(in, "%255s", arquivo) != 1) break;
        }
        // ler coordenadas
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            fscanf(in, "%d %lf %lf", &tempCity, &x[i], &y[i]);
        }

        double *latitude = (double *) malloc((size_t)(N+1) * sizeof(double));
        double *longitude = (double *) malloc((size_t)(N+1) * sizeof(double));

        CalcLatLong ( x, y, N, latitude, longitude );

        // Calcular Matriz Distancia
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = CalcDistGeo ( latitude, longitude, i, j );
            }
        }

    }

    else if ( strcmp(ewt, "ATT") == 0 ) {

        while ( strcmp(arquivo, "NODE_COORD_SECTION") != 0 ) {
            if (fscanf(in, "%255s", arquivo) != 1) break;
        }

        // ler coordenadas
        int tempCity;
        int *tempX = (int *) malloc((size_t)(N+1) * sizeof(int));
        int *tempY = (int *) malloc((size_t)(N+1) * sizeof(int));

        for ( int i = 1; i < N+1; i++ ) {
            fscanf(in, "%d %d %d", &tempCity, &tempX[i], &tempY[i]);
            x[i]=tempX[i];
            y[i]=tempY[i];
        }

        // Calcular Matriz Distancia (Pesudo-Euclidiana)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = CalcDistAtt ( x, y, i, j );
            }
        }

    }

    else if ( strcmp(ewt, "XRAY1") == 0 ) {
        printf("XRAY1 - Nao suportado!\n"); }

    else if ( strcmp(ewt, "XRAY2") == 0 ) {
        printf("XRAY2 - Nao suportado!\n"); }

    else if ( strcmp(ewt, "SPECIAL") == 0 ) {
        printf("SPECIAL - Nao suportado!\n"); }

    fclose(in);

    *Dimension = N;
    *Mdist = dist;
}

double CalcDistEuc ( double *X, double *Y, int I, int J )
{
    return
    sqrt ( pow ( X[I] - X[J], 2 ) + pow ( Y[I] - Y[J], 2 ) );
}

double CalcDistAtt ( double *X, double *Y, int I, int J )
{
    // Calcula Pseudo Distancia Euclidiana
    double rij, tij, dij;

    rij = sqrt ( ( pow ( X[I] - X[J], 2 ) + pow ( Y[I] - Y[J], 2 ) ) / 10 );
    tij = floor ( rij + 0.5 );

    if ( tij < rij )
        dij = tij + 1;
    else
        dij = tij;

    return dij;
}

void CalcLatLong ( double *X, double *Y, int n, double *latit, double* longit )
{
    double PI = 3.141592, min;
    int deg;

    for ( int i = 1; i < n+1; i++ ) {
        deg = (int) X[i];
        min = X[i] - deg;
        latit[i] = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
    }

    for ( int i = 1; i < n+1; i++ ) {
        deg = (int) Y[i];
        min = Y[i] - deg;
        longit[i] = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
    }
}

double CalcDistGeo ( double *latit, double *longit, int I, int J )
{
    double q1, q2, q3, RRR = 6378.388;

    q1 = cos( longit[I] - longit[J] );
    q2 = cos( latit[I] - latit[J] );
    q3 = cos( latit[I] + latit[J] );

    return
    (int) ( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0);
}
