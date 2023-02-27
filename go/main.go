package main

import (
    //"bufio"
    "fmt"
    "log"
    "math"
	"math/rand"
    "os"
)

type tSubseq struct {
    T   float64
    W   float64
    C   float64
}

type tInfo struct {
    c [][]float64
    dimen int
    rnd []int
    rnd_index int
}

type tSolution struct {
    s   []int
    seq [][]tSubseq
    cost float64
}

const (
    SWAP          = 0
    REINSERTION    = 1
    OR_OPT_2      = 2
    OR_OPT_3      = 3
    TWO_OPT       = 4
)

func read_data() (int, [][]float64, []int) {
    var (
        cost [][]float64
        dimen   int
        c       int
        rnd_size int
    )

    file, err := os.Open("../distance_matrix")
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()


    fmt.Fscanf(file, "%d", &dimen)

    fmt.Println(dimen)

    cost = make([][]float64, dimen)

    for i := 0; i < dimen; i++ {
        cost[i] = make([]float64, dimen)
        for j := i+1; j < dimen; j++ {
            fmt.Fscanf(file, "%d", &c)

            cost[i][j] = float64(c)
            fmt.Printf("%f ", cost[i][j])

        }
        // o ultimo valor da linha é lido duas vezes para que a função inicie a leitura da próxima linha do arquivo.
        fmt.Fscanf(file, "%d", &c)
        fmt.Println()
    }

    for i := 0; i < dimen; i++ {
        for j := i+1; j < dimen; j++ {
            cost[j][i] = cost[i][j]
        }
    }

    var buff string
    fmt.Fscanf(file, "%s\n", &buff)
    fmt.Fscanf(file, "%s\n", &buff)
    fmt.Fscanf(file, "%s\n", &buff)
    fmt.Println(buff)

    /*
    scanner := bufio.NewScanner(file)
    scanner.Scan()
    scanner.Scan()
    scanner.Scan()
    fmt.Println(scanner.Text())
    */
    //fmt.Fscanf(file, "%s", &buff)
    //fmt.Fscanf(file, "%s", &buff)


    fmt.Fscanf(file, "%d", &rnd_size)
    fmt.Println(rnd_size)

    rnd  := make([]int, rnd_size)

    for i := 0; i < rnd_size; i++ {
        fmt.Fscanf(file, "%d", &rnd[i])
        //if i < 10 {fmt.Println(rnd[i])}
    }

    return dimen, cost, rnd

}

func remove (arr []int, i int) []int {
    arr[i] =  arr[len(arr)-1]
    return arr[:len(arr)-1]
}

func sort(arr *[]int, r int, info tInfo) {

	for i := 0; i < len(*arr); i++ {
        for j := 0; j < len(*arr)-i-1; j++ {
            if info.c[r][(*arr)[j]] > info.c[r][(*arr)[j+1]] {
                tmp := (*arr)[j]
                (*arr)[j] = (*arr)[j+1]
                (*arr)[j+1] = tmp
            }
        }
    }
}

func construction(alpha float64, info *tInfo) []int {
    s := make([]int, 1)
	s[0] = 0

    cL := make([]int, info.dimen-1)
    for i := 0; i < info.dimen-1; i++ {
        cL[i] = i+1
    }

	r := 0
    for len(cL) > 0 {
		sort(&cL, r, *info)

		rg := int(math.Ceil(float64(len(cL)) * alpha))

		index := rand.Intn(rg)
        index = info.rnd[info.rnd_index];
        info.rnd_index++

		c := cL[index]
		r = c
		
		cL = remove(cL, index)
		s = append(s, c)

		//fmt.Println(cL)
		//fmt.Println(s)
    }

	s = append(s, 0)

    return s
}

func NewSolution(info tInfo) tSolution {
	solut := tSolution {}

	solut.s = make([]int, info.dimen+1)
	solut.seq = make([][]tSubseq, info.dimen+1)
    for i := 0; i < info.dimen+1; i++ {
        solut.seq[i] = make([]tSubseq, info.dimen+1)
    }

    return solut

}

func subseq_load(solut *tSolution, info tInfo) {

    for i := 0; i < info.dimen+1; i++ {
        k := 1 - i;

        if i == 0 {
            k--
        }

        solut.seq[i][i].T = 0.0
        solut.seq[i][i].C = 0.0
        if i == 0 {
            solut.seq[i][i].W = 0.0
        } else {
            solut.seq[i][i].W = 1.0
        }

        for j := i+1; j < info.dimen+1; j++ {
            j_prev := j-1
            
            T := info.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev].T
            solut.seq[i][j].T = T

            C := solut.seq[i][j].T + solut.seq[i][j_prev].C
            solut.seq[i][j].C = C

            W := float64(j + k)
            solut.seq[i][j].W = W

        }
    }

    solut.cost = solut.seq[0][info.dimen].C

}

func swap(solut *tSolution, i int, j int) {
    tmp := solut.s[i]
    solut.s[i] = solut.s[j]
    solut.s[j] = tmp
}

func search_swap(solut *tSolution, info tInfo) bool {
    
    var cost_concat_1  float64
    var cost_concat_2  float64
    var cost_concat_3  float64
    var cost_concat_4  float64

    var cost_best float64 = math.MaxFloat64
    var cost_new  float64
    I := 0
    J := 0

    for i :=  1; i < info.dimen-1; i++ {
        i_prev := i - 1
        i_next := i + 1


        cost_concat_1 =   solut.seq[0][ i_prev].T + info.c[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[i][ i_next].T + info.c[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[0][ i_prev].C +
            solut.seq[i][ i_next].W     * cost_concat_1 + info.c[solut.s[i_next]][solut.s[i]] +
            solut.seq[i_next+1][ info.dimen].W * cost_concat_2 + solut.seq[i_next+1][info.dimen ].C

        if cost_new < cost_best {
            cost_best = cost_new
            I = i
            J = i_next
        }

        for j :=  i_next+1; j < info.dimen; j++ {
            j_next := j+1
            j_prev := j-1

            cost_concat_1 = solut.seq[0][i_prev].T + info.c[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + info.c[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][ j_prev].T + info.c[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3  + info.c[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[0][ i_prev].C +
                    cost_concat_1 +
                    solut.seq[i_next][j_prev].W * cost_concat_2 + solut.seq[i_next][ j_prev].C +
                    cost_concat_3 +
                    solut.seq[j_next][info.dimen].W * cost_concat_4 + solut.seq[j_next][info.dimen].C 


            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }


    if cost_best < solut.seq[0][info.dimen].C {
        //println!("swap \n{}", cost_best);
        swap(solut, I, J)
        subseq_load(solut, info)
        //subseq_load(s, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true
    } else {
        return false
    }
}

func reverse(solut * tSolution, i int, j int) {
    f := i
    l := j

    m := int((i+j) /2)

    for f < m {
        swap(solut, f, l)
        f++
        l--
    }
}

func search_two_opt(solut  *tSolution, info tInfo) bool {
    var cost_new  float64
    cost_best := math.MaxFloat64

    var cost_concat_1 float64
    var cost_concat_2 float64

    I := 0
    J := 0

    for i := 1; i < info.dimen-1; i++ {
        i_prev := i - 1;
        rev_seq_cost := solut.seq[i][i+1].T

        for j := i+2; j < info.dimen; j++ {
            j_next := j + 1

            rev_seq_cost += info.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][ j].W-1.0)

            cost_concat_1 =  solut.seq[0][ i_prev].T + info.c[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][ j].T + info.c[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[0][i_prev].C +
                    solut.seq[i][j].W      * cost_concat_1 + rev_seq_cost +
                    solut.seq[j_next][ info.dimen].W * cost_concat_2 + solut.seq[j_next][ info.dimen].C;

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }


    if cost_best < solut.cost {
        reverse(solut, I, J)
        subseq_load(solut, info)
        return true
    } else {
        return false
    }
}

func reinsert(solut * tSolution, i int, j int, pos int) {
    sub := make([]int, j-i+1)
    copy(sub, solut.s[i:j+1])
    if pos < i {

    } else {
    }

}


func search_reinsertion(solut * tSolution, info tInfo, opt int) bool {
    cost_best := f64::MAX
    var cost_new float64

    var cost_concat_1 float64
    var cost_concat_2 float64
    var cost_concat_3 float64

    I := 0;
    J := 0;
    POS := 0;

    for i := 1; i < info.dimen-opt+1; i++ {
        j := opt+i-1
        i_prev := i-1
        j_next := j+1

        for k := 0 k < i_prev; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][k].T + info.c[solut.s[k]][solut.s[i]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][j].T + info.c[solut.s[j]][solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev].T + info.c[solut.s[i_prev]][solut.s[j_next]];

              cost_new = solut.seq[0][k].C +                                                            /*        1st subseq */
                solut.seq[i][j].W              * cost_concat_1 + solut.seq[i][j].C  +                 /* concat 2nd subseq (reinserted seq) */
                solut.seq[k_next][i_prev].W   * cost_concat_2 + solut.seq[k_next][ i_prev].C  +       /* concat 3rd subseq */
                solut.seq[j_next][ info.dimen].W * cost_concat_3 + solut.seq[j_next][ info.dimen].C    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
                POS = k
            }
        }

        for k := i+opt; k < info.dimen; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][ i_prev].T + info.c[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[j_next][ k].T + info.c[solut.s[k]][solut.s[i]]
            cost_concat_3 = cost_concat_2 + solut.seq[i][ j].T + info.c[solut.s[j]][solut.s[k_next]]

            cost_new = solut.seq[0][ i_prev].C  +                                                       /*      1st subseq */
                solut.seq[j_next][k].W         * cost_concat_1 + solut.seq[j_next][ k].C  +           /* concat 2nd subseq */
                solut.seq[i][ j].W              * cost_concat_2 + solut.seq[i][ j].C   +              /* concat 3rd subseq (reinserted seq) */
                solut.seq[k_next][ info.dimen].W * cost_concat_3 + solut.seq[k_next][ info.dimen].C    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
                POS = k
            }
        }
    }


    if cost_best < solut.cost {
        reinsert(solut, I, J, POS+1)
        subseq_load(solut, info)
        return true
    } else {
        return false
    }
}

func RVND(info tInfo) {
    n_list_b := []int{SWAP, REINSERTION, OR_OPT_2, OR_OPT_3, TWO_OPT}
    n_list := make([]int, 5)

    copy(n_list, n_list_b)
    improve := false

    for len(n_list) > 0 {
        index := rand.Intn(len(n_list))
        index = info.rnd[info.rnd_index]
        info.rnd_index++

        switch n_list[index] {
        case REINSERTION:
        case SWAP:
        case OR_OPT_2:
        case OR_OPT_3:
        case TWO_OPT:
        }

        if improve == true {
            copy(n_list, n_list_b)

        } else {
            n_list = remove(n_list, index)
        }
    }
}

func GILS_RVND(Imax int, Iils int , R []float64, info tInfo) {

	solut_crnt := NewSolution(info)
	solut_partial := NewSolution(info)
	solut_best := NewSolution(info)

    for i := 0; i < Imax; i++ {
        index := rand.Intn(len(R))
        index = info.rnd[info.rnd_index]
        info.rnd_index++

        solut_crnt.s = construction(R[index], &info)
        subseq_load(&solut_crnt, info)

        solut_partial.s = solut_crnt.s
        solut_partial.cost = solut_crnt.cost

        iterILS := 0
        for iterILS < Iils {
            //RVND
            if solut_crnt.cost < solut_partial.cost {
                solut_partial.s = solut_crnt.s
                solut_partial.cost = solut_crnt.cost
            }

            // perturbación
            iterILS++
        }

        if solut_partial.cost < solut_best.cost {
            solut_partial.cost = solut_best.cost
            solut_partial.s = solut_best.s
        }

    }
}

func main() {
    fmt.Println("Hello Vourld!")

    info := tInfo {
        rnd_index: 0,
    }

    info.dimen, info.c, info.rnd = read_data()

	solut := NewSolution(info)

    R := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    index := info.rnd[info.rnd_index]
    info.rnd_index++

	solut.s = construction(R[index], &info)
    fmt.Println(solut.s)

    subseq_load(&solut, info)


    fmt.Println(solut.cost)

}
