package main

import (
    "fmt"
    "math/rand"
    "os"
)

type tSeqInfo struct {
    T, W, C float64
}

type tRnd struct {
    rnd []int
    index int
}

var dimension int
var cost [][]float64

func ternary(s bool, t int, f int) int {if s {return t} else {return f}}

func subseq_load(s []int, seq [][]tSeqInfo) {

    for i := 0; i < dimension+1; i++ {
        k := 1 - i - ternary(i == 0, 1, 0)

        seq[i][i].T = 0.0
        seq[i][i].C = 0.0
        seq[i][i].W = float64(ternary(i != 0, 1, 0))
        for j := i+1; j < dimension+1; j++ {
            var j_prev = j-1
            seq[i][j].T = cost[s[j_prev]][s[j]] + seq[i][j_prev].T
            seq[i][j].C = seq[i][j].T + seq[i][j_prev].C
            seq[i][j].W = float64(j + k)
        }
    }
}

func remove(s []int, i int) {
    return append(s[:i], s[i+1:]...)
}

func swap(s []int, i int, j int) {
    var tmp = s[i]
    s[i] = s[j]
    s[j] = tmp
}

func sort(arr []int, r int) {
    for i := 0; i < len(arr); i++ {
        for j := 0; j < len(arr)-i-1; j++ {
            if cost[r][arr[j]] > cost[r][arr[j+1]] {
                swap(arr, j, j+1)
            }
        }
    }
}

func construct(alpha float64, rnd * tRnd) []int {
    s := make([]int, 0)

    cL := make([]int, dimension-1)
    for i := 1; i < dimension; i++ {
        cL[i-1] = i
    }

    var r = 0
    for len(cL) > 0 {
        sort(cL, r)
        fmt.Println(r)

        var rang =  int(float64(len(cL)) * alpha + 1.0)
        var index = rand.Intn(rang)
        //r_index = rnd.index; rnd.index++
        //index = rnd.rnd[r_index]
        var c = cL[index]
        r = c
        fmt.Println(cL)
        cL = remove(cL, index)
        fmt.Println(cL)
        s = append(s, c)
        os.Exit(0)
    }

    fmt.Println(s)

    return s
}

//func GILS_RVND(Imax int, Iils int, R []float64) {
func GILS_RVND(rnd tRnd) {
    Imax := 10
    Iils := ternary(dimension < 100, dimension, 100)
    _ = Iils
    R := []float64{0.0, 0.1, 0.2}

    var s_best []int
    var s_crnt []int
    var s_partial []int

    _, _, _ = s_best, s_crnt, s_partial

    seq := make([][]tSeqInfo, dimension+1)
    for i := 0; i < dimension+1; i++ {
        seq[i] = make([]tSeqInfo, dimension+1)
    }

    var alpha = R[rand.Intn(len(R))]
    construct(alpha, &rnd)
    for i := 0; i < Imax; i++ {

    }
}

func main() {
    var rnd tRnd
    dimension, cost, rnd.rnd = loadData()
    _ = dimension
    _ = cost
    _ = rnd

    fmt.Println(dimension)
    fmt.Println(cost)

    GILS_RVND(rnd)
}
