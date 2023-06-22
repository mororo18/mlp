package main

import (
    //"bufio"
    "time"
    "fmt"
    "log"
    "math"
	"math/rand"
    "os"
)

/*
type tSubseq struct {
    T   float64
    W   float64
    C   float64
}
*/

const (
    T = 0
    W = 1
    C = 2
)

type tInfo struct {
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
    //seq [][]tSubseq
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

func feasible_(s []int, info tInfo) bool {
    is := make([]bool, info.dimen)

    for i:=0; i < info.dimen; i++ {
        is[i] = false
    }

    for i:=0; i < info.dimen; i++ {
        is[s[i]] = true
    }

    for i:=0; i < info.dimen; i++ {
        if is[i] == false {
            return false
        }
    }

    //fmt.Println(s)

    return true
}

func feasible(solut * tSolution, info tInfo) bool {
    is := make([]bool, info.dimen)

    for i:=0; i < info.dimen; i++ {
        is[i] = false
    }

    for i:=0; i < info.dimen; i++ {
        is[solut.s[i]] = true
    }

    for i:=0; i < info.dimen; i++ {
        if is[i] == false {
            return false
        }
    }

    //fmt.Println(solut.s)

    return true
}

func calc_cost(solut * tSolution, info tInfo) float64 {
    total := 0.0 
    n := info.dimen

    for i := 0; i < info.dimen; i++ {
        total += info.c[solut.s[i]][solut.s[i+1]] * float64(n)
        n--
    }

    return total
}

func remove (arr []int, i int) []int {
    return append(arr[:i], arr[i+1:]...)
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
	solut.seq = make([][][]float64, info.dimen+1)
    for i := 0; i < info.dimen+1; i++ {
        solut.seq[i] = make([][]float64, info.dimen+1)
        for j := 0; j < info.dimen+1; j++ {
            solut.seq[i][j] = make([]float64, 3)
        }
    }

    return solut

}

func subseq_load(solut *tSolution, info tInfo) {

    for i := 0; i < info.dimen+1; i++ {
        k := 1 - i;

        if i == 0 {
            k--
        }

        solut.seq[i][i][info.T] = 0.0
        solut.seq[i][i][info.C] = 0.0
        if i == 0 {
            solut.seq[i][i][info.W] = 0.0
        } else {
            solut.seq[i][i][info.W] = 1.0
        }

        for j := i+1; j < info.dimen+1; j++ {
            j_prev := j-1
            
            T := info.c[solut.s[j_prev]][solut.s[j]] +
            solut.seq[i][j_prev][info.T]
            solut.seq[i][j][info.T] = T

            C := solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C]
            solut.seq[i][j][info.C] = C

            W := float64(j + k)
            solut.seq[i][j][info.W] = W

        }
    }

    solut.cost = solut.seq[0][info.dimen][info.C]

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


        cost_concat_1 =   solut.seq[0][ i_prev][info.T] + info.c[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[i][ i_next][info.T] + info.c[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[0][ i_prev][info.C] +
            solut.seq[i][ i_next][info.W]     * cost_concat_1 + info.c[solut.s[i_next]][solut.s[i]] +
            solut.seq[i_next+1][ info.dimen][info.W] * cost_concat_2 + solut.seq[i_next+1][info.dimen ][info.C]

        if cost_new < cost_best {
            cost_best = cost_new
            I = i
            J = i_next
        }

        for j :=  i_next+1; j < info.dimen; j++ {
            j_next := j+1
            j_prev := j-1

            cost_concat_1 = solut.seq[0][i_prev][info.T] + info.c[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + info.c[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][ j_prev][info.T] + info.c[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3  + info.c[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[0][ i_prev][info.C] +
                    cost_concat_1 +
                    solut.seq[i_next][j_prev][info.W] * cost_concat_2 + solut.seq[i_next][ j_prev][info.C] +
                    cost_concat_3 +
                    solut.seq[j_next][info.dimen][info.W] * cost_concat_4 + solut.seq[j_next][info.dimen][info.C] 




            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }

    if cost_best < solut.seq[0][info.dimen][info.C] {
        //println!("swap \n{}", cost_best);
        swap(solut, I, J)

        if (feasible(solut, info) == false) {
            fmt.Println("qebro swap\n")
            os.Exit(0)
        }

        subseq_load(solut, info)

        //fmt.Println(calc_cost(solut, info), cost_best)
        if (calc_cost(solut, info) != cost_best) {
            fmt.Println("qebro swap\n")
            os.Exit(0)
        }

        //fmt.Println("swap", solut.cost)
        //subseq_load(s, info);
        //println!("{}", seq[0][info.dimension][C]);
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

func search_two_opt(solut  *tSolution, info tInfo) bool {
    var cost_new  float64
    cost_best := math.MaxFloat64

    var cost_concat_1 float64
    var cost_concat_2 float64

    I := 0
    J := 0

    for i := 1; i < info.dimen-1; i++ {
        i_prev := i - 1;
        rev_seq_cost := solut.seq[i][i+1][info.T]

        for j := i+2; j < info.dimen; j++ {
            j_next := j + 1

            rev_seq_cost += info.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][ j][info.W]-1.0)

            cost_concat_1 =  solut.seq[0][ i_prev][info.T] + info.c[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][ j][info.T] + info.c[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[0][i_prev][info.C] +
                    solut.seq[i][j][info.W]      * cost_concat_1 + rev_seq_cost +
                    solut.seq[j_next][ info.dimen][info.W] * cost_concat_2 + solut.seq[j_next][ info.dimen][info.C];

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
            }
        }
    }


    if cost_best < solut.cost {
        //fmt.Println(solut.s)
        //fmt.Println(solut.s, I, J)
        reverse(solut, I, J)

        //antes := solut.cost

        if (feasible(solut, info) == false) {
            fmt.Println("qebro two_opt")
            os.Exit(0)
        }

        subseq_load(solut, info)

        if (calc_cost(solut, info) != cost_best) {
            //fmt.Println(solut.s, "qebro two_opt")
            //fmt.Println("Antes ", antes, "\nDepois ", solut.cost)
            fmt.Println("Cost Best ", cost_best)
            os.Exit(0)
        }


        //fmt.Println("two_opt", solut.cost)

        return true
    } 

    return false
}

func reinsert(solut * tSolution, i int, j int, pos int) {
    sz := j-i+1
    sub := make([]int, sz)

    /*
    remove_seq := func(arr []int, a int, b int) []int {
        return append(arr[:a], arr[b+1:]...)
    }
    */
     //s_cpy := make([]int, len(solut.s))

    copy(sub, solut.s[i:j+1])
    //copy(s_cpy, solut.s)

    if pos < i {
        //fmt.Println("Antes", solut.s)
        copy(solut.s[pos+sz:j+1], solut.s[pos:i])
        copy(solut.s[pos:pos+sz], sub)
        //solut.s = remove_seq(solut.s, i, j)
        //sub = append(solut.s[:pos], sub...)
        //solut.s = append(sub, solut.s[pos:]...)

    } else {
        //fmt.Println("Depois", solut.s)
        copy(solut.s[i:i+pos-j], solut.s[j+1:pos])
        copy(solut.s[pos-(j-i+1): pos], sub)
        //sub = append((*solut).s[:pos], sub...)
        //fmt.Println("sub=",sub, "solut.s", solut.s)
        //sub = append(sub, solut.s[pos+1:]...)
        //fmt.Println("sub=",sub)
        //solut.s = remove_seq(solut.s, i, j)
        //fmt.Println("solut.s=",solut.s)
    }

}


func search_reinsertion(solut * tSolution, info tInfo, opt int) bool {
    cost_best := math.MaxFloat64
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

        for k := 0; k < i_prev; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][k][info.T] + info.c[solut.s[k]][solut.s[i]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T] + info.c[solut.s[j]][solut.s[k_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][info.T] + info.c[solut.s[i_prev]][solut.s[j_next]];

              cost_new = solut.seq[0][k][info.C] +                                                            /*        1st subseq */
                solut.seq[i][j][info.W]              * cost_concat_1 + solut.seq[i][j][info.C]  +                 /* concat 2nd subseq (reinserted seq) */
                solut.seq[k_next][i_prev][info.W]   * cost_concat_2 + solut.seq[k_next][ i_prev][info.C]  +       /* concat 3rd subseq */
                solut.seq[j_next][ info.dimen][info.W] * cost_concat_3 + solut.seq[j_next][ info.dimen][info.C]    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new
                I = i
                J = j
                POS = k
            }
        }

        for k := i+opt; k < info.dimen; k++ {
            k_next := k+1

            cost_concat_1 = solut.seq[0][ i_prev][info.T] + info.c[solut.s[i_prev]][solut.s[j_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[j_next][ k][info.T] + info.c[solut.s[k]][solut.s[i]]
            cost_concat_3 = cost_concat_2 + solut.seq[i][ j][info.T] + info.c[solut.s[j]][solut.s[k_next]]

            cost_new = solut.seq[0][ i_prev][info.C]  +                                                       /*      1st subseq */
                solut.seq[j_next][k][info.W]         * cost_concat_1 + solut.seq[j_next][ k][info.C]  +           /* concat 2nd subseq */
                solut.seq[i][ j][info.W]              * cost_concat_2 + solut.seq[i][ j][info.C]   +              /* concat 3rd subseq (reinserted seq) */
                solut.seq[k_next][ info.dimen][info.W] * cost_concat_3 + solut.seq[k_next][ info.dimen][info.C]    /* concat 4th subseq */

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

        if (feasible(solut, info) == false) {
            //fmt.Println("qebro reinsert\n")
            os.Exit(0)
        }

        if (calc_cost(solut, info) != cost_best) {
            fmt.Println("qebro reinsert\n", solut.s)
            os.Exit(0)
        }

        subseq_load(solut, info)
        //fmt.Println("reinsert", solut.cost)

        return true
    }

    return false
}

func RVND(solut * tSolution, info * tInfo) {
    n_list_b := []int{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
    n_list := make([]int, 5)

    copy(n_list, n_list_b)

    //fmt.Println(solut.s)

    for len(n_list) > 0 {
        index := rand.Intn(len(n_list))
        index = info.rnd[info.rnd_index]
        info.rnd_index++

        //fmt.Println(n_list)

        improve := false
        switch n_list[index] {
        case REINSERTION:
            improve = search_reinsertion(solut, *info, REINSERTION)
        case OR_OPT_2:
            improve = search_reinsertion(solut, *info, OR_OPT_2)
        case OR_OPT_3:
            improve = search_reinsertion(solut, *info, OR_OPT_3)
        case TWO_OPT:
            improve = search_two_opt(solut, *info)
        case SWAP:
            improve = search_swap(solut, *info)
        }

        if improve == true {
            n_list = make([]int, 5)
            copy(n_list, n_list_b)

        } else {
            n_list = remove(n_list, index)
        }
    }
}

func perturb(sl []int, info * tInfo) []int {
    //fmt.Println("Perturbacion")
    s := make([]int, info.dimen+1)
    copy(s, sl)

    A_start := 1
    A_end   := 1
    B_start := 1
    B_end   := 1

    size_max := 3
    if int(float64(len(s) / 10.0)) >= 2 {
        size_max = int(float64(len(s) / 10.0))
    }

    size_min := 2


    for (A_start <= B_start &&  B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {
        A_start = rand.Intn(len(s) - 2 - size_max) + 1
        A_end = A_start + rand.Intn(size_max - size_min) + size_min

        B_start = rand.Intn(len(s) - 2 - size_max) + 1
        B_end = B_start + rand.Intn(size_max - size_min) + size_min

        A_start = info.rnd[info.rnd_index]
        info.rnd_index++
        A_end = A_start + info.rnd[info.rnd_index]
        info.rnd_index++

        B_start = info.rnd[info.rnd_index]
        info.rnd_index++
        B_end = B_start + info.rnd[info.rnd_index]
        info.rnd_index++

    }

    reinsert_s := func (s_ * []int, i int, j int, pos int) {
        sz := j-i+1
        sub := make([]int, j-i+1)

        copy(sub, (*s_)[i:j+1])
        //fmt.Println("sub", sub)

        if pos < i {
            //fmt.Println("Antes ", s_)
            copy((*s_)[pos+sz:j+1], (*s_)[pos:i])
            //fmt.Println(s_)
            copy((*s_)[pos:pos+sz], sub)
            //fmt.Println(s_)
            //solut.s = remove_seq(solut.s, i, j)
            //sub = append(solut.s[:pos], sub...)
            //solut.s = append(sub, solut.s[pos:]...)

        } else {
            //fmt.Println("Depois", s_)
            copy((*s_)[i:i+pos-j], (*s_)[j+1:pos])
            copy((*s_)[pos-(j-i+1): pos], sub)
        }

        if feasible_(*s_, *info) == false {
            //fmt.Println("Perturb qebrad")
            os.Exit(0)
        }


    }

    //fmt.Println("B s e", B_start, B_end)



    if A_start < B_start {
        reinsert_s(&s, B_start, B_end-1, A_end);
        //fmt.Println("Reinsercao 1", s)
        reinsert_s(&s, A_start, A_end-1, B_end);
        //fmt.Println("Reinsercao 2", s)
    } else {
        reinsert_s(&s, A_start, A_end-1, B_end);
        //fmt.Println("Reinsercao 1", s)
        reinsert_s(&s, B_start, B_end-1, A_end);
        //fmt.Println("Reinsercao 2", s)
    }

    //fmt.Println("PERTURBACOA Resultado = ", s)

    return s;
}

func GILS_RVND(Imax int, Iils int , R [26]float64, info tInfo) {

	solut_crnt := NewSolution(info)
	solut_partial := NewSolution(info)
	solut_best := NewSolution(info)
    solut_best.cost = math.Inf(1)

    for i := 0; i < Imax; i++ {

        fmt.Printf("[+] Local Search %d\n", i)

        index := rand.Intn(len(R))
        index = info.rnd[info.rnd_index]
        info.rnd_index++

        solut_crnt.s = construction(R[index], &info)
        subseq_load(&solut_crnt, info)
        fmt.Printf("\t[+] Constructing Inital Solution.. %.2f\n", solut_crnt.cost)
        fmt.Println("\t", solut_crnt.s)

        copy(solut_partial.s, solut_crnt.s)
        solut_partial.cost = solut_crnt.cost

        fmt.Println("\t[+] Looking for the best Neighbor..")
        iterILS := 0
        for iterILS < Iils {
            RVND(&solut_crnt, &info)
            if solut_crnt.cost < solut_partial.cost {
                copy(solut_partial.s, solut_crnt.s)
                solut_partial.cost = solut_crnt.cost
                iterILS = 0
            }


            solut_crnt.s = perturb(solut_partial.s, &info)
            subseq_load(&solut_crnt, info)
            //fmt.Println(solut_crnt.cost, solut_crnt.s)
            // perturbación
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
    fmt.Println("Hello Vourld!")

    info := tInfo {
        rnd_index: 0,
    }

    info.dimen, info.c, info.rnd = read_data()
    info.T = T
    info.W = W
    info.C = C

	//solut := NewSolution(info)

    R := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    Imax := 10
    Iils := 100
    if info.dimen < 100 {
        Iils = info.dimen
    }

    start := time.Now()

    GILS_RVND(Imax, Iils, R, info)

    t := time.Now()
    elapsed := t.Sub(start)
    fmt.Println("TIME:", elapsed.Seconds())
    //GILS_RVND()

}
