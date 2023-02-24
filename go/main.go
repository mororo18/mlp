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
    cost [][]float64
    dimen int
    rnd []int
    rnd_index int
}

type tSolution struct {
    s   []int
    seq [][]tSubseq
    cost float64
}

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

func sort(arr *[]int, r int, info tInfo) {

	for i := 0; i < len(*arr); i++ {
        for j := 0; j < len(*arr)-i-1; j++ {
            if info.cost[r][(*arr)[j]] > info.cost[r][(*arr)[j+1]] {
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

    remove := func(arr []int, i int) []int {
        arr[i] =  arr[len(arr)-1]
        return arr[:len(arr)-1]
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
            
            T := info.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev].T
            solut.seq[i][j].T = T

            C := solut.seq[i][j].T + solut.seq[i][j_prev].C;
            solut.seq[i][j].C = C

            W := float64(j + k)
            solut.seq[i][j].W = W

        }
    }

    solut.cost = solut.seq[0][info.dimen].C;

}

func GILS_RVND(Imax int, Iils int , R []float64) {

	solut_crnt := NewSolution(info)
	solut_partial := NewSolution(info)
	solut_best := NewSolution(info)

    for i := 0; i < Imax; i++ {
        index := info.rnd[info.rnd_index]
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

    info.dimen, info.cost, info.rnd = read_data()

	solut := NewSolution(info)

    R := [...]float64{0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    index := info.rnd[info.rnd_index]
    info.rnd_index++

	solut.s = construction(R[index], &info)
    fmt.Println(solut.s)

    subseq_load(&solut, info)


    fmt.Println(solut.cost)

}
