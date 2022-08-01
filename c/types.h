#ifndef TYPES
#define TYPES

#include <string.h>
#include <float.h>
#include <stdlib.h>

typedef struct tInfo {
    double ** cost;
    int dimen;
    int * rnd;
    int rnd_index;
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
static int to_1D(const int i, const int j, const int size);
#endif

/*==========================SET==========================*/

static void seq_set_C(tSolution * solut, int i, int j, double value);
static void seq_set_T(tSolution * solut, int i, int j, double value);
static void seq_set_W(tSolution * solut, int i, int j, double value);

/*==========================GET==========================*/
static double seq_get_C(const tSolution * solut, int i, int j);
static double seq_get_T(const tSolution * solut, int i, int j);
static double seq_get_W(const tSolution * solut, int i, int j);
tSolution Solution_init(tInfo info);
void      Solution_cpy(tSolution * src, tSolution * tgt, const tInfo * info);
void      Solution_free(tSolution * solut);

void tInfo_free(tInfo * info);

/*==========================inline==========================*/

#ifdef FLAT
inline
int to_1D(const int i, const int j, const int size) {
    return i * size + j;
}
#endif


inline
void seq_set_C(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].C = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].C = value;
#endif
}

inline
void seq_set_T(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].T = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].T = value;
#endif
}

inline
void seq_set_W(tSolution * solut, int i, int j, double value) {
#ifdef MATRIX
    solut->seq[i][j].W = value;
#elif defined(FLAT)
    solut->seq[to_1D(i, j, solut->s_size)].W = value;
#endif
}


inline
double seq_get_C(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].C;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].C;
#endif
}

inline
double seq_get_T(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].T;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].T;
#endif
}

inline 
double seq_get_W(const tSolution * solut, int i, int j) {
#ifdef MATRIX
    return solut->seq[i][j].W;
#elif defined(FLAT)
    return solut->seq[to_1D(i, j, solut->s_size)].W;
#endif
}

#endif
