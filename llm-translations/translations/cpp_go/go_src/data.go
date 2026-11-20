package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

func matrixAllocate(sz int) [][]float64 {
	ptr := make([][]float64, sz)

	for i := 0; i < sz; i++ {
		ptr[i] = make([]float64, sz)
	}
	return ptr
}

func loadData(matrix *[][]float64, rnd *[]int) int {
	var dimension int

	file, err := os.Open("../distance_matrix")

	if err == nil {
		reader := bufio.NewReader(file)

		line, _ := readLine(reader)
		dimension, _ = strconv.Atoi(line)
		*matrix = matrixAllocate(dimension)

		flag := false
		rndSize := -1
		i := 0
		for {
			line, ok := readLine(reader)
			if !ok {
				break
			}

			j := i + 1
			for strings.Contains(line, " ") {
				pos := strings.Index(line, " ")
				value, _ := strconv.ParseFloat(line[0:pos], 64)
				line = line[pos+1:]
				(*matrix)[i][j] = value
				(*matrix)[j][i] = value
				j++
			}

			i++

			if flag {

				if rndSize < 0 {
					rndSize, _ = strconv.Atoi(line)
					*rnd = append(make([]int, 0, rndSize), (*rnd)...)
				} else {
					v, _ := strconv.Atoi(line)
					*rnd = append(*rnd, v)
				}
			}

			if line == "RND" {
				flag = true
			}
		}
	}

	if file != nil {
		file.Close()
	}
	return dimension
}

// readLine mimics std::getline: returns the line without its trailing newline
// and false once there is nothing left to read.
func readLine(reader *bufio.Reader) (string, bool) {
	line, err := reader.ReadString('\n')
	if err != nil && line == "" {
		return "", false
	}
	line = strings.TrimRight(line, "\n")
	line = strings.TrimRight(line, "\r")
	return line, true
}

type Data struct {
	dimen int
	cost  [][]float64
	rnd   []int
}

func (d *Data) load() {
	d.dimen = loadData(&d.cost, &d.rnd)
}

func (d *Data) getRndVec() []int {
	return d.rnd
}

func (d *Data) getDimen() int {
	return d.dimen
}

func (d *Data) getCostPtr() [][]float64 {
	return d.cost
}
