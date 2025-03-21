package main

import (
    "time"
    "fmt"
    "log"
    "math"
	"math/rand"
    "os"
)

const (
    T = 0
    W = 1
    C = 2
)

type tData struct {
    T int
    W int
    C int
    c [][]float64
    dimen int
    rnd []int
    rnd_index int
}

type tSolution struct {
    s   []int
    seq [][][]float64
    cost float64
}

const (
    SWAP          = 0
    REINSERTION   = 1
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

        }
        // o ultimo valor da linha é lido duas vezes para que a função inicie a leitura da próxima linha do arquivo.
        fmt.Fscanf(file, "%d", &c)
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

    fmt.Fscanf(file, "%d", &rnd_size)

    rnd  := make([]int, rnd_size)

    for i := 0; i < rnd_size; i++ {
        fmt.Fscanf(file, "%d", &rnd[i])
    }

    return dimen, cost, rnd
}

func remove (arr []int, i int) []int {
    return append(arr[:i], arr[i+1:]...)
}

func sort(arr *[]int, r int, data *tData) {
	quicksort(arr, 0, len(*arr)-1, data, r)
}

func quicksort(arr *[]int, left, right int, data *tData, r int) {
	if left < right {
		pivot := partition(arr, left, right, data, r)
		quicksort(arr, left, pivot-1, data, r)
		quicksort(arr, pivot+1, right, data, r)
	}
}

func partition(arr *[]int, left, right int, data *tData, r int) int {
	pivot := (*arr)[right]
	i := left - 1
	for j := left; j < right; j++ {
		if data.c[r][(*arr)[j]] < data.c[r][pivot] {
			i++
			(*arr)[i], (*arr)[j] = (*arr)[j], (*arr)[i]
		}
	}
	(*arr)[i+1], (*arr)[right] = (*arr)[right], (*arr)[i+1]
	return i + 1
}

func construction(alpha float64, data *tData) []int {
    s := make([]int, 1)
	s[0] = 0

    cL := make([]int, data.dimen-1)
    for i := 0; i < data.dimen-1; i++ {
        cL[i] = i+1
    }

	r := 0
    for len(cL) > 0 {
		sort(&cL, r, data)

        index := data.rnd[data.rnd_index];
        data.rnd_index++

		c := cL[index]
		r = c
		
		cL = remove(cL, index)
		s = append(s, c)
    }

	s = append(s, 0)

    return s
}

func NewSolution(data tData) tSolution {
	solut := tSolution {}

	solut.s = make([]int, data.dimen+1)
	solut.seq = make([][][]float64, data.dimen+1)
    for i := 0; i < data.dimen+1; i++ {
        solut.seq[i] = make([][]float64, data.dimen+1)
        for j := 0; j < data.dimen+1; j++ {
            solut.seq[i][j] = make([]float64, 3)
        }
    }

    return solut

}

func update_subseq_info_matrix(solut *tSolution, data *tData) {

    for i := 0; i < data.dimen+1; i++ {
        k := 1 - i;

        if i == 0 {
            k--
        }

        solut.seq[i][i][data.T] = 0.0
        solut.seq[i][i][data.C] = 0.0
        if i == 0 {
            solut.seq[i][i][data.W] = 0.0
        } else {
            solut.seq[i][i][data.W] = 1.0
        }

        for j := i+1; j < data.dimen+1; j++ {
            j_prev := j-1
            
            T := data.c[solut.s[j_prev]][solut.s[j]] +
            solut.seq[i][j_prev][data.T]
            solut.seq[i][j][data.T] = T

            C := solut.seq[i][j][data.T] + solut.seq[i][j_prev][data.C]
            solut.seq[i][j][data.C] = C

            W := float64(j + k)
            solut.seq[i][j][data.W] = W

        }
    }

    solut.cost = solut.seq[0][data.dimen][data.C]

}

func swap(solut *tSolution, i int, j int) {
    tmp := solut.s[i]
    solut.s[i] = solut.s[j]
    solut.s[j] = tmp
}

func search_swap(solut *tSolution, data *tData) bool {
    
    var cost_concat_1  float64
    var cost_concat_2  float64
    var cost_concat_3  float64
    var cost_concat_4  float64

    var cost_best float64 = math.MaxFloat64
    var cost_new  float64
    I := 0
    J := 0

    for i :=  1; i < data.dimen-1; i++ {
        i_prev := i - 1
        i_next := i + 1


        cost_concat_1 =   solut.seq[0][ i_prev][data.T] + data.c[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[i][ i_next][data.T] + data.c[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[0][ i_prev][data.C] +
            solut.seq[i][ i_next][data.W]     * cost_concat_1 + data.c[solut.s[i_next]][solut.s[i]] +
            solut.seq[i_next+1][ data.dimen][data.W] * cost_concat_2 + solut.seq[i_next+1][data.dimen ][data.C]

        if cost_new < cost_best {
            cost_best = cost_new
            I = i
            J = i_next
        }

        for j :=  i_next+1; j < data.dimen; j++ {
            j_next := j+1
            j_prev := j-1

            cost_concat_1 = solut.seq[0][i_prev][data.T] + data.c[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + data.c[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][ j_prev][data.T] + data.c[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3  + data.c[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[0][ i_prev][data.C] +
                    cost_concat_1 +
                    solut.seq[i_next][j_prev][data.W] * cost_concat_2 + solut.seq[i_next][ j_prev][data.C] +
                    cost_concat_3 +
                    solut.seq[j_next][data.dimen][data.W] * cost_concat_4 + solut.seq[j_next][data.dimen][data.C] 

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }

    if cost_best < solut.seq[0][data.dimen][data.C] {
        swap(solut, I, J)
        update_subseq_info_matrix(solut, data)
        return true
    }

    return false
}

func reverse(solut * tSolution, i int, j int) {
    f := i
    l := j

    m := int((i+j) /2)

    for f <= m {
        swap(solut, f, l)
        f++
        l--
    }
}

func search_two_opt(solut  *tSolution, data *tData) bool {
    var cost_new  float64
    cost_best := math.MaxFloat64

    var cost_concat_1 float64
    var cost_concat_2 float64

    I := 0
    J := 0

    for i := 1; i < data.dimen-1; i++ {
        i_prev := i - 1;
        rev_seq_cost := solut.seq[i][i+1][data.T]

        for j := i+2; j < data.dimen; j++ {
            j_next := j + 1

            rev_seq_cost += data.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][ j][data.W]-1.0)

            cost_concat_1 =  solut.seq[0][ i_prev][data.T] + data.c[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][ j][data.T] + data.c[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[0][i_prev][data.C] +
                    solut.seq[i][j][data.W]      * cost_concat_1 + rev_seq_cost +
                    solut.seq[j_next][ data.dimen][data.W] * cost_concat_2 + solut.seq[j_next][ data.dimen][data.C];

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }


    if cost_best < solut.cost {
        reverse(solut, I, J)
        update_subseq_info_matrix(solut, data)
        return true
    } 

    return false
}

func reinsert(solut * tSolution, i int, j int, pos int) {
    sz := j-i+1
    sub := make([]int, sz)

    copy(sub, solut.s[i:j+1])

    if pos < i {
        copy(solut.s[pos+sz:j+1], solut.s[pos:i])
        copy(solut.s[pos:pos+sz], sub)

    } else {
        copy(solut.s[i:i+pos-j], solut.s[j+1:pos])
        copy(solut.s[pos-(j-i+1): pos], sub)
    }

}


func search_reinsertion(solut * tSolution, data *tData, opt int) bool {
    cost_best := math.MaxFloat64
    var cost_new float64

    var cost_concat_1 float64
    var cost_concat_2 float64
    var cost_concat_3 float64

    I := 0;
    J := 0;
    POS := 0;

    for i := 1; i < data.dimen-opt+1; i++ {
        j := opt+i-1
        i_prev := i-1
        j_next := j+1

        for k := 0; k < i_prev; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][k][data.T] + data.c[solut.s[k]][solut.s[i]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][j][data.T] + data.c[solut.s[j]][solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][data.T] + data.c[solut.s[i_prev]][solut.s[j_next]];

              cost_new = solut.seq[0][k][data.C] +                                                            /*        1st subseq */
                solut.seq[i][j][data.W]              * cost_concat_1 + solut.seq[i][j][data.C]  +                 /* concat 2nd subseq (reinserted seq) */
                solut.seq[k_next][i_prev][data.W]   * cost_concat_2 + solut.seq[k_next][ i_prev][data.C]  +       /* concat 3rd subseq */
                solut.seq[j_next][ data.dimen][data.W] * cost_concat_3 + solut.seq[j_next][ data.dimen][data.C]    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
                POS = k
            }
        }

        for k := i+opt; k < data.dimen; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][ i_prev][data.T] + data.c[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[j_next][ k][data.T] + data.c[solut.s[k]][solut.s[i]]
            cost_concat_3 = cost_concat_2 + solut.seq[i][ j][data.T] + data.c[solut.s[j]][solut.s[k_next]]

            cost_new = solut.seq[0][ i_prev][data.C]  +                                                       /*      1st subseq */
                solut.seq[j_next][k][data.W]         * cost_concat_1 + solut.seq[j_next][ k][data.C]  +           /* concat 2nd subseq */
                solut.seq[i][ j][data.W]              * cost_concat_2 + solut.seq[i][ j][data.C]   +              /* concat 3rd subseq (reinserted seq) */
                solut.seq[k_next][ data.dimen][data.W] * cost_concat_3 + solut.seq[k_next][ data.dimen][data.C]    /* concat 4th subseq */

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
        update_subseq_info_matrix(solut, data)

        return true
    }

    return false
}

func RVND(solut * tSolution, data * tData) {
    n_list_b := []int{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
    n_list := make([]int, 5)

    copy(n_list, n_list_b)

    for len(n_list) > 0 {
        index := data.rnd[data.rnd_index]
        data.rnd_index++

        improve := false
        switch n_list[index] {
        case REINSERTION:
            improve = search_reinsertion(solut, data, REINSERTION)
        case OR_OPT_2:
            improve = search_reinsertion(solut, data, OR_OPT_2)
        case OR_OPT_3:
            improve = search_reinsertion(solut, data, OR_OPT_3)
        case TWO_OPT:
            improve = search_two_opt(solut, data)
        case SWAP:
            improve = search_swap(solut, data)
        }

        if improve == true {
            n_list = make([]int, 5)
            copy(n_list, n_list_b)

        } else {
            n_list = remove(n_list, index)
        }
    }
}

func perturb(sl []int, data * tData) []int {
    s := make([]int, data.dimen+1)
    copy(s, sl)

    A_start := 1
    A_end   := 1
    B_start := 1
    B_end   := 1

    for (A_start <= B_start &&  B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {

        A_start = data.rnd[data.rnd_index]
        data.rnd_index++
        A_end = A_start + data.rnd[data.rnd_index]
        data.rnd_index++

        B_start = data.rnd[data.rnd_index]
        data.rnd_index++
        B_end = B_start + data.rnd[data.rnd_index]
        data.rnd_index++

    }

    reinsert_s := func (s_ * []int, i int, j int, pos int) {
        sz := j-i+1
        sub := make([]int, j-i+1)

        copy(sub, (*s_)[i:j+1])

        if pos < i {
            copy((*s_)[pos+sz:j+1], (*s_)[pos:i])
            copy((*s_)[pos:pos+sz], sub)

        } else {
            copy((*s_)[i:i+pos-j], (*s_)[j+1:pos])
            copy((*s_)[pos-(j-i+1): pos], sub)
        }

    }

    if A_start < B_start {
        reinsert_s(&s, B_start, B_end-1, A_end);
        reinsert_s(&s, A_start, A_end-1, B_end);
    } else {
        reinsert_s(&s, A_start, A_end-1, B_end);
        reinsert_s(&s, B_start, B_end-1, A_end);
    }

    return s;
}

func GILS_RVND(Imax int, Iils int , R [26]float64, data tData) {

	solut_crnt := NewSolution(data)
	solut_partial := NewSolution(data)
	solut_best := NewSolution(data)
    solut_best.cost = math.Inf(1)

    for i := 0; i < Imax; i++ {

        fmt.Printf("[+] Local Search %d\n", i)

        index := rand.Intn(len(R))
        index = data.rnd[data.rnd_index]
        data.rnd_index++

        solut_crnt.s = construction(R[index], &data)
        update_subseq_info_matrix(&solut_crnt, &data)
        fmt.Printf("\t[+] Constructing Inital Solution.. %.2f\n", solut_crnt.cost)
        fmt.Println("\t", solut_crnt.s)

        copy(solut_partial.s, solut_crnt.s)
        solut_partial.cost = solut_crnt.cost

        fmt.Println("\t[+] Looking for the best Neighbor..")
        iterILS := 0
        for iterILS < Iils {
            RVND(&solut_crnt, &data)
            if solut_crnt.cost < solut_partial.cost {
                copy(solut_partial.s, solut_crnt.s)
                solut_partial.cost = solut_crnt.cost
                iterILS = 0
            }


            solut_crnt.s = perturb(solut_partial.s, &data)
            update_subseq_info_matrix(&solut_crnt, &data)
            iterILS++
        }

        if solut_partial.cost < solut_best.cost {
            solut_best.cost = solut_partial.cost
            copy(solut_best.s, solut_partial.s)
        }

    }

    fmt.Println("COST: ", solut_best.cost)
    fmt.Println("SOLUTION: ", solut_best.s)
}

func main() {
    data := tData {
        rnd_index: 0,
    }

    data.dimen, data.c, data.rnd = read_data()
    data.T = T
    data.W = W
    data.C = C

    R := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    Imax := 10
    Iils := 100
    if data.dimen < 100 {
        Iils = data.dimen
    }

    start := time.Now()

    GILS_RVND(Imax, Iils, R, data)

    t := time.Now()
    elapsed := t.Sub(start)
    fmt.Println("TIME:", elapsed.Seconds())
}
