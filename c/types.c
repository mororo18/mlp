#include "types.h"

/*==========================SET==========================*/


/*==========================GET==========================*/

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

void Solution_free(tSolution * solut) {
    free(solut->s);

#ifdef MATRIX
    for (int i = 0; i < solut->s_size; i++) {
        free(solut->seq[i]);
    }
#endif
    free(solut->seq);
}

void tInfo_free(tInfo * info) {
    for (int i = 0; i < info->dimen; i++) {
        free(info->cost[i]);
    }
    free(info->cost);
    free(info->rnd);
}

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

