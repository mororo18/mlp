package main

import (
    "fmt"
    "os"
    "log"
    "bufio"
    "strings"
)

func loadData() {
    var dimension int
    var line string
    var value float64
    var cost [][]float64

    fmt.Println("Hello, World!")

    f, err := os.Open("../distance_matrix")

    if err != nil {
        log.Fatal(err)
    }

    defer f.Close()

    scanner := bufio.NewScanner(f)

    scanner.Scan(); line = scanner.Text()
    fmt.Sscanf(line, "%d ", &dimension)
    fmt.Println(dimension)

    cost = make([][]float64, dimension)

    for i:=0; i<dimension; i++ {
        cost[i] = make([]float64, dimension)
    }


    /*
    scanner.Scan(); line = scanner.Text()

    fmt.Println(line)
    fmt.Println(line[0:2])
    fmt.Println(line[0])
    fmt.Println()
    */

    for i := 0; i < dimension; i++ {
        scanner.Scan(); line = scanner.Text()
        strs := strings.Split(line, " ")
        index := len(strs)-1
        strs = append(strs[:index], strs[index+1:]...)
        //fmt.Printf("%q\n", strs);
        for j := i+1; j < dimension; j++ {
            //index :=  strings.Index(line, " ")
            fmt.Sscanf(strs[j-(i+1)], "%v", &value)
            cost[i][j] = value;
            cost[j][i] = value;
        }
        fmt.Println( cost[i])
    }

    /*
    scanner.Scan(); line = scanner.Text()
    fmt.Println(line)

    fmt.Sscanf(line, "%v ", &value)
    fmt.Println(value)
    for scanner.Scan() {
        fmt.Println(scanner.Text())
    }
    */

    if err :=  scanner.Err(); err != nil {
        log.Fatal(err)
    }
}
