package main

import (
    "fmt"
)

var dimension int
var cost [][]float64
var rnd []int

type tSeqInfo struct {
    T, W, C float64
}

func subseq_load(s []int, seq [][]tSeqInfo) {
    ternary := func(s bool, t int, f int) int {if s {return t} else {return f}}

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

func main() {
    dimension, cost, rnd = loadData()
    _ = dimension
    _ = cost
    _ = rnd

    fmt.Println(dimension)
    fmt.Println(cost)


    s := make([]int, dimension+1)
    for i := 0; i < dimension; i++ {
        s[i] = i
    }
    seq := make([][]tSeqInfo, dimension+1)
    for i := 0; i < dimension+1; i++ {
        seq[i] = make([]tSeqInfo, dimension+1)
    }
    s[dimension] = 0
    fmt.Println(s)

    subseq_load(s, seq)
    fmt.Println(seq[0][dimension].C)
}
