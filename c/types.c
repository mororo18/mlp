#ifndef TYPES
#define TYPES

#include <string.h>
#include <float.h>
#include <stdlib.h>

typedef unsigned uint;

typedef struct tInfo {
    double ** cost;
    int dimen;
    uint T;
    uint C;
    uint W;
    int * rnd;
    uint rnd_index;
} tInfo;

typedef struct tSeqInfo {
    double T, C, W;
} tSeqInfo;

typedef tSeqInfo tSeq;
typedef tSeq * tSeq_;
typedef tSeq_ * tSeq__;

typedef struct tSolution {
#ifdef MATRIX
    tSeq__ seq;
#elif defined(FLAT)
    tSeq_ seq;
#endif
    //double *** seq;
    double cost;
    int * s;
    int s_size;
    size_t size;
    float MBsize;
    //int s_size;
} tSolution;

#ifdef FLAT
static
inline
int to_1D(const int i, const int j, const int size) {
    return i * size + j;
}
#endif

/*==========================SET==========================*/

static
inline
void seq_set_C(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].C = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].C = value;
#endif
}

static
inline
void seq_set_T(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].T = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].T = value;
#endif
}

static
inline
void seq_set_W(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].W = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].W = value;
#endif
}


/*==========================GET==========================*/
static
inline
double seq_get_C(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].C;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].C;
#endif
}

static
inline
double seq_get_T(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].T;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].T;
#endif
}

static
inline 
double seq_get_W(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].W;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].W;
#endif
}

static
tSolution Solution_init(tInfo info) {
    tSolution solut;
    solut.size = sizeof(tSolution);
    solut.s = (int*) calloc(info.dimen+1, sizeof(int));
    //solut.s_size = info.dimen+1;

  //solut.seq = (double ***) calloc(info.dimen+1, sizeof(double **));
  //for (int i = 0; i < info.dimen+1; i++) {
  //    solut.seq[i] = (double **) calloc(info.dimen+1, sizeof(double *));
  //    for (int j = 0; j < info.dimen+1; j++) {
  //        solut.seq[i][j] = (double *) calloc(3, sizeof(double));
  //    }
  //}

    solut.s_size = info.dimen + 1;

#ifdef MATRIX
    solut.seq = (tSeq__) calloc(info.dimen+1, sizeof(tSeq_));
    for (int i = 0; i < info.dimen+1; i++) {
        solut.seq[i] = (tSeq_) calloc(info.dimen+1, sizeof(tSeq));
    }

    solut.size += solut.s_size * sizeof(tSeq_);
    solut.size += solut.s_size * solut.s_size * sizeof(tSeq);
#elif defined(FLAT)
    solut.seq = (tSeq_) calloc((info.dimen+1)*(info.dimen+1), sizeof(tSeq));

    solut.size += solut.s_size * solut.s_size * sizeof(tSeq);
#endif

    solut.MBsize = solut.size / (1024.0 * 1024);

    solut.cost = DBL_MAX;

    return solut;
}

static
void Solution_free(tSolution * solut) {
    free(solut->s);

#ifdef MATRIX
    for (int i = 0; i < solut->s_size; i++) {
        free(solut->seq[i]);
    }
#endif
    free(solut->seq);
}

static
void tInfo_free(tInfo * info) {
    for (int i = 0; i < info->dimen; i++) {
        free(info->cost[i]);
    }
    free(info->cost);
    free(info->rnd);
}

static
void Solution_cpy(tSolution * src, tSolution * tgt, const tInfo * info) {

    memcpy(tgt->s, src->s, sizeof(int)*(info->dimen+1));
    tgt->cost = src->cost;

    /*
    for (int i = 0; i < info.dimen+1; i++) {
        for (int j = 0; j < info.dimen+1; j++) {
            //memcpy(tgt.seq[i][j], src.seq[i][j], 3 * sizeof(double));
            std::copy(src.seq[i][j], src.seq[i][j] + 3, tgt.seq[i][j]);
        }
    }
    */

}


#endif
