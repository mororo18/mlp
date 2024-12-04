#include "types.h"

/*==========================SET==========================*/


/*==========================GET==========================*/

tSolution Solution_init(tData data) {
    tSolution solut;
    solut.size = sizeof(tSolution);
    solut.s = (int*) calloc(data.dimen+1, sizeof(int));

    solut.s_size = data.dimen + 1;

#ifdef MATRIX
    solut.seq = (tSeq__) calloc(data.dimen+1, sizeof(tSeq_));
    for (int i = 0; i < data.dimen+1; i++) {
        solut.seq[i] = (tSeq_) calloc(data.dimen+1, sizeof(tSeq));
    }

    solut.size += solut.s_size * sizeof(tSeq_);
    solut.size += solut.s_size * solut.s_size * sizeof(tSeq);
#elif defined(FLAT)
    solut.seq = (tSeq_) calloc((data.dimen+1)*(data.dimen+1), sizeof(tSeq));

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

void tData_free(tData * data) {
    for (int i = 0; i < data->dimen; i++) {
        free(data->cost[i]);
    }
    free(data->cost);
    free(data->rnd);
}

void Solution_cpy(tSolution * src, tSolution * tgt, const tData * data) {

    memcpy(tgt->s, src->s, sizeof(int)*(data->dimen+1));
    tgt->cost = src->cost;

    /*
    for (int i = 0; i < info.dimen+1; i++) {
        for (int j = 0; j < info.dimen+1; j++) {
            std::copy(src.seq[i][j], src.seq[i][j] + 3, tgt.seq[i][j]);
        }
    }
    */

}

