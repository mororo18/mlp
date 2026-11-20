package main

import (
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"syscall"
	"time"
)

const (
	TABLE_SZ = 26
	DBL_SZ   = 8
	INT_SZ   = 4
)

const dblEpsilon = 2.220446049250313e-16
const dblMax = math.MaxFloat64

type tData struct {
	cost      [][]float64
	dimen     int
	rnd       []int
	rnd_index uint
}

type tSeqInfo struct {
	T float64
	C float64
	W float64
}

type Neighborhood int

const (
	REINSERTION Neighborhood = iota + 1
	OR_OPT_2
	OR_OPT_3
	SWAP
	TWO_OPT
)

type tSolution struct {
	s    []int
	seq  []tSeqInfo
	cost float64
}

func to_1D(x, y, size int) int {
	a := x*(size+1) + y

	return a
}

func Solution_init(data tData) tSolution {
	var solut tSolution
	solut.s = make([]int, data.dimen+1)

	solut.seq = make([]tSeqInfo, (data.dimen+1)*(data.dimen+1))
	solut.cost = dblMax

	return solut
}

func Solution_cpy(src *tSolution, tgt *tSolution, data *tData) {

	tgt.s = append(tgt.s[:0], src.s...)
	tgt.cost = src.cost

}

func R_table(i int) float64 {
	table := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
		0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}
	return table[i]
}

func print_s(s []int) {

	for i := 0; i < len(s); i++ {
		fmt.Print(s[i]+1, " ")
	}
	fmt.Println()
}

func partition(arr []int, left, right int, data *tData, r int) int {
	pivot := arr[right]
	i := left - 1
	for j := left; j < right; j++ {
		if data.cost[r][arr[j]] < data.cost[r][pivot] {
			i++
			arr[i], arr[j] = arr[j], arr[i]
		}
	}
	arr[i+1], arr[right] = arr[right], arr[i+1]
	return i + 1
}

func quicksort(arr []int, left, right int, data *tData, r int) {
	if left < right {
		pivot := partition(arr, left, right, data, r)
		quicksort(arr, left, pivot-1, data, r)
		quicksort(arr, pivot+1, right, data, r)
	}
}

func sortCL(arr []int, r int, data *tData) {
	quicksort(arr, 0, len(arr)-1, data, r)
}

func construct(alpha float64, data *tData) []int {

	s := make([]int, 0, data.dimen+1)
	s = append(s, 0)
	cL := make([]int, data.dimen-1)
	for i := 1; i < data.dimen; i++ {
		cL[i-1] = i
	}

	r := 0
	for len(cL) != 0 {
		sortCL(cL, r, data)

		index := int(data.rnd[data.rnd_index])
		data.rnd_index++

		c := cL[index]
		s = append(s, c)
		r = c
		cL = append(cL[:index], cL[index+1:]...)
	}

	s = append(s, 0)

	return s
}

func swap(vec []int, i, j int) {
	vec[i], vec[j] = vec[j], vec[i]
}

func reverse(vec []int, i, j int) {
	for i < j {
		vec[i], vec[j] = vec[j], vec[i]
		i++
		j--
	}
}

func reinsert(vec []int, i, j, pos int) {
	seq := make([]int, j-i+1)
	copy(seq, vec[i:j+1])
	if pos < i {
		copy(vec[pos+len(seq):], vec[pos:i])
		copy(vec[pos:], seq)
	} else {
		copy(vec[i:], vec[j+1:pos])
		copy(vec[pos-len(seq):], seq)
	}

}

func update_subseq_info_matrix(solut *tSolution, data *tData, index int) {
	var i, j, j_prev, k int
	from := index
	var t bool
	for i = 0; i < data.dimen+1; i++ {
		k = 1 - i
		if i == 0 {
			k = 1 - i - 1
		}
		t = i == from

		solut.seq[to_1D(i, i, data.dimen)].T = 0.0
		solut.seq[to_1D(i, i, data.dimen)].C = 0.0
		if i == 0 {
			solut.seq[to_1D(i, i, data.dimen)].W = 0.0
		} else {
			solut.seq[to_1D(i, i, data.dimen)].W = 1.0
		}

		for j = i + 1; j < data.dimen+1; j++ {
			j_prev = j - 1

			solut.seq[to_1D(i, j, data.dimen)].T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[to_1D(i, j_prev, data.dimen)].T
			solut.seq[to_1D(i, j, data.dimen)].C = solut.seq[to_1D(i, j, data.dimen)].T + solut.seq[to_1D(i, j_prev, data.dimen)].C
			solut.seq[to_1D(i, j, data.dimen)].W = float64(j + k)
		}
		if t {
			from++
		}
	}

	solut.cost = solut.seq[to_1D(0, data.dimen, data.dimen)].C
}

func search_swap(solut *tSolution, data *tData) bool {
	var cost_new, cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4 float64
	cost_best := dblMax
	var i, j, j_prev, j_next, i_prev, i_next int
	var I int
	var J int

	for i = 1; i < data.dimen-1; i++ {
		i_prev = i - 1
		i_next = i + 1

		//consecutive nodes
		cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[i_next]]
		cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, data.dimen)].T + data.cost[solut.s[i]][solut.s[i_next+1]]

		cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C + //       1st subseq
			solut.seq[to_1D(i, i_next, data.dimen)].W*(cost_concat_1) + data.cost[solut.s[i_next]][solut.s[i]] + // concat 2nd subseq
			solut.seq[to_1D(i_next+1, data.dimen, data.dimen)].W*(cost_concat_2) + solut.seq[to_1D(i_next+1, data.dimen, data.dimen)].C // concat 3rd subseq

		if cost_new < cost_best {
			cost_best = cost_new - dblEpsilon
			I = i
			J = i_next
		}

		for j = i_next + 1; j < data.dimen; j++ {
			j_next = j + 1
			j_prev = j - 1

			cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j]]
			cost_concat_2 = cost_concat_1 + data.cost[solut.s[j]][solut.s[i_next]]
			cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, data.dimen)].T + data.cost[solut.s[j_prev]][solut.s[i]]
			cost_concat_4 = cost_concat_3 + data.cost[solut.s[i]][solut.s[j_next]]

			cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C + // 1st subseq
				cost_concat_1 + // concat 2nd subseq (single node)
				solut.seq[to_1D(i_next, j_prev, data.dimen)].W*cost_concat_2 + solut.seq[to_1D(i_next, j_prev, data.dimen)].C + // concat 3rd subseq
				cost_concat_3 + // concat 4th subseq (single node)
				solut.seq[to_1D(j_next, data.dimen, data.dimen)].W*cost_concat_4 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C // concat 5th subseq

			if cost_new < cost_best {
				cost_best = cost_new - dblEpsilon
				I = i
				J = j
			}
		}
	}

	if cost_best < solut.cost-dblEpsilon {
		swap(solut.s, I, J)
		update_subseq_info_matrix(solut, data, I)
		return true
	}

	return false
}

func search_two_opt(solut *tSolution, data *tData) bool {
	var cost_new, cost_concat_1, cost_concat_2 float64
	cost_best := dblMax
	var rev_seq_cost float64
	var i, j, i_prev, j_next int
	var I int
	var J int

	for i = 1; i < data.dimen-1; i++ {
		i_prev = i - 1

		rev_seq_cost = solut.seq[to_1D(i, i+1, data.dimen)].T
		for j = i + 2; j < data.dimen; j++ {
			j_next = j + 1

			rev_seq_cost += data.cost[solut.s[j-1]][solut.s[j]] * (solut.seq[to_1D(i, j, data.dimen)].W - 1.0)

			cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[j]][solut.s[i_prev]]
			cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j_next]][solut.s[i]]

			cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C + //  1st subseq
				solut.seq[to_1D(i, j, data.dimen)].W*cost_concat_1 + rev_seq_cost + // concat 2nd subseq (reversed seq)
				solut.seq[to_1D(j_next, data.dimen, data.dimen)].W*cost_concat_2 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C // concat 3rd subseq

			if cost_new < cost_best {
				cost_best = cost_new - dblEpsilon
				I = i
				J = j
			}
		}
	}

	if cost_best < solut.cost-dblEpsilon {
		reverse(solut.s, I, J)
		update_subseq_info_matrix(solut, data, 0)
		return true
	}

	return false
}

func search_reinsertion(solut *tSolution, data *tData, opt int) bool {
	var cost_new, cost_concat_1, cost_concat_2, cost_concat_3 float64
	cost_best := dblMax
	var i, j, k, k_next, i_prev, j_next int
	var I int
	var J int
	var POS int

	for i, j = 1, opt+1-1; i < data.dimen-opt+1; i, j = i+1, j+1 {
		j_next = j + 1
		i_prev = i - 1

		//k -> edges
		for k = 0; k < i_prev; k++ {
			k_next = k + 1

			cost_concat_1 = solut.seq[to_1D(0, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]]
			cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]
			cost_concat_3 = cost_concat_2 + solut.seq[to_1D(k_next, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]

			cost_new = solut.seq[to_1D(0, k, data.dimen)].C + //       1st subseq
				solut.seq[to_1D(i, j, data.dimen)].W*cost_concat_1 + solut.seq[to_1D(i, j, data.dimen)].C + //  concat 2nd subseq (reinserted seq)
				solut.seq[to_1D(k_next, i_prev, data.dimen)].W*cost_concat_2 + solut.seq[to_1D(k_next, i_prev, data.dimen)].C + //  concat 3rd subseq
				solut.seq[to_1D(j_next, data.dimen, data.dimen)].W*cost_concat_3 + solut.seq[to_1D(j_next, data.dimen, data.dimen)].C // concat 4th subseq

			if cost_new < cost_best {
				cost_best = cost_new - dblEpsilon
				I = i
				J = j
				POS = k
			}
		}

		for k = i + opt; k < data.dimen; k++ {
			k_next = k + 1

			cost_concat_1 = solut.seq[to_1D(0, i_prev, data.dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]
			cost_concat_2 = cost_concat_1 + solut.seq[to_1D(j_next, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]]
			cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]

			cost_new = solut.seq[to_1D(0, i_prev, data.dimen)].C + //       1st subseq
				solut.seq[to_1D(j_next, k, data.dimen)].W*cost_concat_1 + solut.seq[to_1D(j_next, k, data.dimen)].C + // concat 2nd subseq
				solut.seq[to_1D(i, j, data.dimen)].W*cost_concat_2 + solut.seq[to_1D(i, j, data.dimen)].C + // concat 3rd subseq (reinserted seq)
				solut.seq[to_1D(k_next, data.dimen, data.dimen)].W*cost_concat_3 + solut.seq[to_1D(k_next, data.dimen, data.dimen)].C // concat 4th subseq

			if cost_new < cost_best {
				cost_best = cost_new - dblEpsilon
				I = i
				J = j
				POS = k
			}
		}
	}

	if cost_best < solut.cost-dblEpsilon {
		reinsert(solut.s, I, J, POS+1)
		from := I
		if POS+1 < I {
			from = POS + 1
		}
		update_subseq_info_matrix(solut, data, from)
		return true
	}

	return false
}

func RVND(solut *tSolution, data *tData) {

	neighbd_list := []Neighborhood{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
	var index uint
	var neighbd Neighborhood
	var improve_flag bool

	for len(neighbd_list) != 0 {

		index = uint(data.rnd[data.rnd_index])
		data.rnd_index++
		neighbd = neighbd_list[index]
		improve_flag = false

		switch neighbd {
		case REINSERTION:
			improve_flag = search_reinsertion(solut, data, int(REINSERTION))
		case OR_OPT_2:
			improve_flag = search_reinsertion(solut, data, int(OR_OPT_2))
		case OR_OPT_3:
			improve_flag = search_reinsertion(solut, data, int(OR_OPT_3))
		case SWAP:
			improve_flag = search_swap(solut, data)
		case TWO_OPT:
			improve_flag = search_two_opt(solut, data)
		}
		if improve_flag {
			neighbd_list = []Neighborhood{SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
		} else {
			neighbd_list = append(neighbd_list[:index], neighbd_list[index+1:]...)
		}

	}

}

func perturb(solut *tSolution, data *tData) []int {
	s := append([]int(nil), solut.s...)
	A_start := 1
	A_end := 1
	B_start := 1
	B_end := 1

	size_max := int(math.Floor(float64(data.dimen+1) / 10))
	if size_max < 2 {
		size_max = 2
	}
	size_min := 2
	_ = size_max
	_ = size_min
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

	if A_start < B_start {
		reinsert(s, B_start, B_end-1, A_end)
		reinsert(s, A_start, A_end-1, B_end)
	} else {
		reinsert(s, A_start, A_end-1, B_end)
		reinsert(s, B_start, B_end-1, A_end)
	}

	return s
}

func GILS_RVND(Imax, Iils int, data *tData, verbose bool) {

	solut_partial := Solution_init(*data)
	solut_crnt := Solution_init(*data)
	solut_best := Solution_init(*data)

	for i := 0; i < Imax; i++ {
		aux := data.rnd[data.rnd_index]
		data.rnd_index++

		alpha := R_table(aux)

		if verbose {
			fmt.Printf("[+] Search %d\n", i+1)
			fmt.Printf("\t[+] Constructing..\n")
		}

		solut_crnt.s = construct(alpha, data)
		update_subseq_info_matrix(&solut_crnt, data, 0)

		Solution_cpy(&solut_crnt, &solut_partial, data)
		if verbose {
			fmt.Printf("\t[+] Looking for the best Neighbor..\n")
			fmt.Printf("\t    Construction Cost: %.3f\n", solut_partial.cost)
		}

		iterILS := 0
		for iterILS < Iils {
			RVND(&solut_crnt, data)
			if solut_crnt.cost < solut_partial.cost-dblEpsilon {
				Solution_cpy(&solut_crnt, &solut_partial, data)
				iterILS = 0
			}

			solut_crnt.s = perturb(&solut_partial, data)
			update_subseq_info_matrix(&solut_crnt, data, 0)
			iterILS++
		}

		if solut_partial.cost < solut_best.cost-dblEpsilon {
			Solution_cpy(&solut_partial, &solut_best, data)
		}

		if verbose {
			fmt.Println("\tCurrent best cost: " + fmtDouble(solut_best.cost))

			fmt.Print("SOLUCAO: ")
			for i := 0; i < len(solut_best.s); i++ {
				fmt.Print(solut_best.s[i], " ")
			}
			fmt.Println()
		}

	}
	fmt.Printf("COST: %.2f\n", solut_best.cost)
}

// fmtDouble reproduces the default formatting of std::ostream for doubles
// (6 significant digits, shortest of fixed/scientific).
func fmtDouble(v float64) string {
	s := strconv.FormatFloat(v, 'G', 6, 64)
	if i := strings.IndexByte(s, 'E'); i >= 0 {
		mant := s[:i]
		exp := s[i+1:]
		sign := "+"
		if exp[0] == '+' || exp[0] == '-' {
			sign = string(exp[0])
			exp = exp[1:]
		}
		for len(exp) < 2 {
			exp = "0" + exp
		}
		s = mant + "e" + sign + exp
	}
	return s
}

func cpuTime() float64 {
	var ru syscall.Rusage
	syscall.Getrusage(syscall.RUSAGE_SELF, &ru)
	user := float64(ru.Utime.Sec) + float64(ru.Utime.Usec)/1e6
	sys := float64(ru.Stime.Sec) + float64(ru.Stime.Usec)/1e6
	return user + sys
}

func main() {
	Imax := 10
	var Iils int
	verbose := false

	for i := 1; i < len(os.Args); i++ {
		if os.Args[i] == "-v" || os.Args[i] == "--verbose" {
			verbose = true
		}
	}

	data := tData{}

	var rnd []int

	data.dimen = loadData(&data.cost, &rnd)
	data.rnd = rnd
	data.rnd_index = 0

	_ = Solution_init(data)

	if data.dimen < 100 {
		Iils = data.dimen
	} else {
		Iils = 100
	}

	t1 := time.Now()
	start := cpuTime()
	GILS_RVND(Imax, Iils, &data, verbose)
	cpu_time := cpuTime() - start
	t2 := time.Now()
	duration := t2.Sub(t1).Milliseconds()

	res := float64(duration) / 10e2
	fmt.Println("TIME: " + fmtDouble(cpu_time))
	fmt.Println("wall clock (s): " + fmtDouble(res))

	fmt.Println("Tamanho RND " + strconv.Itoa(len(rnd)))

	os.Exit(0)
}
