package main

import (
    "bufio"
    "fmt"
    "log"
    "os"
)

func main() {
    fmt.Println("Hello Vourld!")

    file, err := os.Open("../distance_matrix")
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()

    var (
        dimen   int
        c    int
    )

    fmt.Fscanf(file, "%d", &dimen)

    fmt.Println(dimen)

    for i := 1; i < dimen; i++ {
        for j := i+1; j < dimen; j++ {
            fmt.Fscanf(file, "%d", &c)

            fmt.Printf("%d ", c)
        }
        fmt.Println()
    }

    scanner := bufio.NewScanner(file)

    if err := scanner.Err(); err != nil {
        log.Fatal(err)
    }
}
