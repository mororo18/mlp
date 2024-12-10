package main

import (
	_ "runtime/pprof"
	"time"
	//"bufio"
	"fmt"
	"log"
	"math"
	"os"
)

type tSubseq struct {
	T float64
	W float64
	C float64
}

type tData struct {
	c         [][]float64
	dimen     int
	rnd       []int
	rnd_index int
}

type tSolution struct {
	s    []int
	seq  [][]tSubseq
	cost float64
}

const (
	SWAP        = 0
	REINSERTION = 1
	OR_OPT_2    = 2
	OR_OPT_3    = 3
	TWO_OPT     = 4
)

func read_data() (int, [][]float64, []int) {
	var (
		cost     [][]float64
		dimen    int
		c        int
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
		for j := i + 1; j < dimen; j++ {
			fmt.Fscanf(file, "%d", &c)

			cost[i][j] = float64(c)
			fmt.Printf("%f ", cost[i][j])

		}
		// o ultimo valor da linha é lido duas vezes para que a função inicie a leitura da próxima linha do arquivo.
		fmt.Fscanf(file, "%d", &c)
		fmt.Println()
	}

	for i := 0; i < dimen; i++ {
		for j := i + 1; j < dimen; j++ {
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

	rnd := make([]int, rnd_size)

	for i := 0; i < rnd_size; i++ {
		fmt.Fscanf(file, "%d", &rnd[i])
		//if i < 10 {fmt.Println(rnd[i])}
	}

	return dimen, cost, rnd

}

func feasible_(s []int, data tData) bool {
	is := make([]bool, data.dimen)

	for i := 0; i < data.dimen; i++ {
		is[i] = false
	}

	for i := 0; i < data.dimen; i++ {
		is[s[i]] = true
	}

	for i := 0; i < data.dimen; i++ {
		if is[i] == false {
			return false
		}
	}

	//fmt.Println(s)

	return true
}

/*
func feasible(solut * tSolution, data tData) bool {
    is := make([]bool, data.dimen)

    for i:=0; i < data.dimen; i++ {
        is[i] = false
    }

    for i:=0; i < data.dimen; i++ {
        is[solut.s[i]] = true
    }

    for i:=0; i < data.dimen; i++ {
        if is[i] == false {
            return false
        }
    }

    //fmt.Println(solut.s)

    return true
}

func calc_cost(solut * tSolution, data tData) float64 {
    total := 0.0
    n := data.dimen

    for i := 0; i < data.dimen; i++ {
        total += data.c[solut.s[i]][solut.s[i+1]] * float64(n)
        n--
    }

    return total
}
*/

func remove(arr []int, i int) []int {
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
		cL[i] = i + 1
	}

	r := 0
	for len(cL) > 0 {
		sort(&cL, r, data)

        index := data.rnd[data.rnd_index]
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
	solut := tSolution{}

	solut.s = make([]int, data.dimen+1)
	solut.seq = make([][]tSubseq, data.dimen+1)
	for i := 0; i < data.dimen+1; i++ {
		solut.seq[i] = make([]tSubseq, data.dimen+1)
	}

	return solut

}

func update_subseq_info_matrix(solut *tSolution, data tData) {

	for i := 0; i < data.dimen+1; i++ {
		k := 1 - i

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

		for j := i + 1; j < data.dimen+1; j++ {
			j_prev := j - 1

			solut.seq[i][j].T = data.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev].T

			solut.seq[i][j].C = solut.seq[i][j].T + solut.seq[i][j_prev].C

			solut.seq[i][j].W = float64(j + k)
			//= W

		}
	}

	solut.cost = solut.seq[0][data.dimen].C

}

func swap(solut *tSolution, i int, j int) {
	tmp := solut.s[i]
	solut.s[i] = solut.s[j]
	solut.s[j] = tmp
}

func search_swap(solut *tSolution, data tData) bool {

	var cost_concat_1 float64
	var cost_concat_2 float64
	var cost_concat_3 float64
	var cost_concat_4 float64

	var cost_best float64 = math.MaxFloat64
	var cost_new float64
	I := 0
	J := 0

	for i := 1; i < data.dimen-1; i++ {
		i_prev := i - 1
		i_next := i + 1

		cost_concat_1 = solut.seq[0][i_prev].T + data.c[solut.s[i_prev]][solut.s[i_next]]
		cost_concat_2 = cost_concat_1 + solut.seq[i][i_next].T + data.c[solut.s[i]][solut.s[i_next+1]]

		cost_new = solut.seq[0][i_prev].C +
			solut.seq[i][i_next].W*cost_concat_1 + data.c[solut.s[i_next]][solut.s[i]] +
			solut.seq[i_next+1][data.dimen].W*cost_concat_2 + solut.seq[i_next+1][data.dimen].C

		if cost_new < cost_best {
			cost_best = cost_new
			I = i
			J = i_next
		}

		for j := i_next + 1; j < data.dimen; j++ {
			j_next := j + 1
			j_prev := j - 1

			cost_concat_1 = solut.seq[0][i_prev].T + data.c[solut.s[i_prev]][solut.s[j]]
			cost_concat_2 = cost_concat_1 + data.c[solut.s[j]][solut.s[i_next]]
			cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev].T + data.c[solut.s[j_prev]][solut.s[i]]
			cost_concat_4 = cost_concat_3 + data.c[solut.s[i]][solut.s[j_next]]

			cost_new = solut.seq[0][i_prev].C +
				cost_concat_1 +
				solut.seq[i_next][j_prev].W*cost_concat_2 + solut.seq[i_next][j_prev].C +
				cost_concat_3 +
				solut.seq[j_next][data.dimen].W*cost_concat_4 + solut.seq[j_next][data.dimen].C

			if cost_new < cost_best {
				cost_best = cost_new
				I = i
				J = j
			}
		}
	}

	if cost_best < solut.seq[0][data.dimen].C {
		//println!("swap \n{}", cost_best);
		swap(solut, I, J)

		update_subseq_info_matrix(solut, data)

		/*
		   if (feasible(solut, data) == false) {
		       fmt.Println("qebro swap\n")
		       os.Exit(0)
		   }

		   //fmt.Println(calc_cost(solut, data), cost_best)
		   if (calc_cost(solut, data) != cost_best) {
		       fmt.Println("qebro swap\n")
		       os.Exit(0)
		   }
		*/

		//fmt.Println("swap", solut.cost)
		//update_subseq_info_matrix(s, data);
		//println!("{}", seq[0][data.dimension][C]);
		return true
	}

	return false
}

func reverse(solut *tSolution, i int, j int) {
	f := i
	l := j

	m := int((i + j) / 2)

	for f <= m {
		swap(solut, f, l)
		f++
		l--
	}
}

func search_two_opt(solut *tSolution, data tData) bool {
	var cost_new float64
	cost_best := math.MaxFloat64

	var cost_concat_1 float64
	var cost_concat_2 float64

	I := 0
	J := 0

	for i := 1; i < data.dimen-1; i++ {
		i_prev := i - 1
		rev_seq_cost := solut.seq[i][i+1].T

		for j := i + 2; j < data.dimen; j++ {
			j_next := j + 1

			rev_seq_cost += data.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][j].W - 1.0)

			cost_concat_1 = solut.seq[0][i_prev].T + data.c[solut.s[j]][solut.s[i_prev]]
			cost_concat_2 = cost_concat_1 + solut.seq[i][j].T + data.c[solut.s[j_next]][solut.s[i]]

			cost_new = solut.seq[0][i_prev].C +
				solut.seq[i][j].W*cost_concat_1 + rev_seq_cost +
				solut.seq[j_next][data.dimen].W*cost_concat_2 + solut.seq[j_next][data.dimen].C

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
		update_subseq_info_matrix(solut, data)

		/*
		   if (feasible(solut, data) == false) {
		       fmt.Println("qebro two_opt")
		       os.Exit(0)
		   }

		   if (calc_cost(solut, data) != cost_best) {
		       //fmt.Println(solut.s, "qebro two_opt")
		       //fmt.Println("Antes ", antes, "\nDepois ", solut.cost)
		       fmt.Println("Cost Best ", cost_best)
		       os.Exit(0)
		   }
		*/

		//fmt.Println("two_opt", solut.cost)

		return true
	}

	return false
}

func reinsert(solut *tSolution, i int, j int, pos int) {
	sz := j - i + 1
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
		copy(solut.s[pos-(j-i+1):pos], sub)
		//sub = append((*solut).s[:pos], sub...)
		//fmt.Println("sub=",sub, "solut.s", solut.s)
		//sub = append(sub, solut.s[pos+1:]...)
		//fmt.Println("sub=",sub)
		//solut.s = remove_seq(solut.s, i, j)
		//fmt.Println("solut.s=",solut.s)
	}

}

func search_reinsertion(solut *tSolution, data tData, opt int) bool {
	cost_best := math.MaxFloat64
	var cost_new float64

	var cost_concat_1 float64
	var cost_concat_2 float64
	var cost_concat_3 float64

	I := 0
	J := 0
	POS := 0

	for i := 1; i < data.dimen-opt+1; i++ {
		j := opt + i - 1
		i_prev := i - 1
		j_next := j + 1

		for k := 0; k < i_prev; k++ {
			k_next := k + 1

			cost_concat_1 = solut.seq[0][k].T + data.c[solut.s[k]][solut.s[i]]
			cost_concat_2 = cost_concat_1 + solut.seq[i][j].T + data.c[solut.s[j]][solut.s[k_next]]
			cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev].T + data.c[solut.s[i_prev]][solut.s[j_next]]

			cost_new = solut.seq[0][k].C + /*        1st subseq */
				solut.seq[i][j].W*cost_concat_1 + solut.seq[i][j].C + /* concat 2nd subseq (reinserted seq) */
				solut.seq[k_next][i_prev].W*cost_concat_2 + solut.seq[k_next][i_prev].C + /* concat 3rd subseq */
				solut.seq[j_next][data.dimen].W*cost_concat_3 + solut.seq[j_next][data.dimen].C /* concat 4th subseq */

			if cost_new < cost_best {
				cost_best = cost_new
				I = i
				J = j
				POS = k
			}
		}

		for k := i + opt; k < data.dimen; k++ {
			k_next := k + 1

			cost_concat_1 = solut.seq[0][i_prev].T + data.c[solut.s[i_prev]][solut.s[j_next]]
			cost_concat_2 = cost_concat_1 + solut.seq[j_next][k].T + data.c[solut.s[k]][solut.s[i]]
			cost_concat_3 = cost_concat_2 + solut.seq[i][j].T + data.c[solut.s[j]][solut.s[k_next]]

			cost_new = solut.seq[0][i_prev].C + /*      1st subseq */
				solut.seq[j_next][k].W*cost_concat_1 + solut.seq[j_next][k].C + /* concat 2nd subseq */
				solut.seq[i][j].W*cost_concat_2 + solut.seq[i][j].C + /* concat 3rd subseq (reinserted seq) */
				solut.seq[k_next][data.dimen].W*cost_concat_3 + solut.seq[k_next][data.dimen].C /* concat 4th subseq */

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

		/*
		   if (feasible(solut, data) == false) {
		       //fmt.Println("qebro reinsert\n")
		       os.Exit(0)
		   }

		   if (calc_cost(solut, data) != cost_best) {
		       fmt.Println("qebro reinsert\n", solut.s)
		       os.Exit(0)
		   }
		*/

		update_subseq_info_matrix(solut, data)
		//fmt.Println("reinsert", solut.cost)

		return true
	}

	return false
}

func RVND(solut *tSolution, data *tData) {
	n_list_b := []int{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
	n_list := make([]int, 5)

	copy(n_list, n_list_b)

	//fmt.Println(solut.s)

	for len(n_list) > 0 {
        index := data.rnd[data.rnd_index]
		data.rnd_index++

		//fmt.Println(n_list)

		improve := false
		switch n_list[index] {
		case REINSERTION:
			improve = search_reinsertion(solut, *data, REINSERTION)
		case OR_OPT_2:
			improve = search_reinsertion(solut, *data, OR_OPT_2)
		case OR_OPT_3:
			improve = search_reinsertion(solut, *data, OR_OPT_3)
		case TWO_OPT:
			improve = search_two_opt(solut, *data)
		case SWAP:
			improve = search_swap(solut, *data)
		}

		if improve == true {
			n_list = make([]int, 5)
			copy(n_list, n_list_b)

		} else {
			n_list = remove(n_list, index)
		}
	}
}

func perturb(sl []int, data *tData) []int {
	//fmt.Println("Perturbacion")
	s := make([]int, data.dimen+1)
	copy(s, sl)

	A_start := 1
	A_end := 1
	B_start := 1
	B_end := 1

	for (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {

		A_start = data.rnd[data.rnd_index]
		data.rnd_index++
		A_end = A_start + data.rnd[data.rnd_index]
		data.rnd_index++

		B_start = data.rnd[data.rnd_index]
		data.rnd_index++
		B_end = B_start + data.rnd[data.rnd_index]
		data.rnd_index++

	}

	reinsert_s := func(s_ *[]int, i int, j int, pos int) {
		sz := j - i + 1
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
			copy((*s_)[pos-(j-i+1):pos], sub)
		}

		if feasible_(*s_, *data) == false {
			//fmt.Println("Perturb qebrad")
			os.Exit(0)
		}

	}

	//fmt.Println("B s e", B_start, B_end)

	if A_start < B_start {
		reinsert_s(&s, B_start, B_end-1, A_end)
		//fmt.Println("Reinsercao 1", s)
		reinsert_s(&s, A_start, A_end-1, B_end)
		//fmt.Println("Reinsercao 2", s)
	} else {
		reinsert_s(&s, A_start, A_end-1, B_end)
		//fmt.Println("Reinsercao 1", s)
		reinsert_s(&s, B_start, B_end-1, A_end)
		//fmt.Println("Reinsercao 2", s)
	}

	//fmt.Println("PERTURBACOA Resultado = ", s)

	return s
}

func GILS_RVND(Imax int, Iils int, R [26]float64, data tData) {

	solut_crnt := NewSolution(data)
	solut_partial := NewSolution(data)
	solut_best := NewSolution(data)
	solut_best.cost = math.Inf(1)

	for i := 0; i < Imax; i++ {

		fmt.Printf("[+] Local Search %d\n", i)

        index := data.rnd[data.rnd_index]
		data.rnd_index++

		solut_crnt.s = construction(R[index], &data)
		update_subseq_info_matrix(&solut_crnt, data)
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
			update_subseq_info_matrix(&solut_crnt, data)
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

	data := tData{
		rnd_index: 0,
	}

	data.dimen, data.c, data.rnd = read_data()

	//solut := NewSolution(data)

	R := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

	Imax := 10
	Iils := 100
	if data.dimen < 100 {
		Iils = data.dimen
	}

	//  f, err := os.Create("privateJourneyFindAll.prof")
	//  if err != nil {
	//      log.Fatal(err)
	//  }
	//  pprof.StartCPUProfile(f)
	//  defer pprof.StopCPUProfile()

	start := time.Now()

	GILS_RVND(Imax, Iils, R, data)

	//elapsed :=
	t := time.Now()
	elapsed := t.Sub(start)
	fmt.Println("TIME:", elapsed.Seconds())
	//GILS_RVND()

}
