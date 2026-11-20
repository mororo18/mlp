package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
)

// tokenReader mimics the token-at-a-time extraction of std::ifstream's
// operator>>.
type tokenReader struct {
	sc *bufio.Scanner
}

func newTokenReader(f *os.File) *tokenReader {
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1024*1024), 1024*1024)
	sc.Split(bufio.ScanWords)
	return &tokenReader{sc: sc}
}

func (t *tokenReader) str() (string, bool) {
	if !t.sc.Scan() {
		return "", false
	}
	return t.sc.Text(), true
}

func (t *tokenReader) int() int {
	s, _ := t.str()
	v, _ := strconv.Atoi(s)
	return v
}

func (t *tokenReader) float() float64 {
	s, _ := t.str()
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

// skipUntil consumes tokens until one of the wanted keywords is read; it
// returns the token that stopped the scan.
func (t *tokenReader) skipUntil(arquivo string, wanted ...string) string {
	for {
		matched := false
		for _, w := range wanted {
			if arquivo == w {
				matched = true
				break
			}
		}
		if matched {
			return arquivo
		}
		s, ok := t.str()
		if !ok {
			return arquivo
		}
		arquivo = s
	}
}

func readData(args []string, Dimension *int, Mdist *[][]float64) {

	if len(args) < 2 {
		fmt.Println("\nFaltando parametros")
		fmt.Println(" ./exec [Instancia] ")
		os.Exit(1)
	}

	if len(args) > 2 {
		fmt.Println("\nMuitos parametros")
		fmt.Println(" ./exec [Instancia] ")
		os.Exit(1)
	}

	var N int
	var arquivo, ewt string

	instancia := args[1]

	in, err := os.Open(instancia)

	if err != nil {
		fmt.Println("arquivo nao pode ser aberto")
		os.Exit(1)
	}
	defer in.Close()

	r := newTokenReader(in)

	arquivo = r.skipUntil(arquivo, "DIMENSION:", "DIMENSION")

	if arquivo == "DIMENSION" {
		arquivo, _ = r.str()
	}

	N = r.int()

	arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_TYPE:", "EDGE_WEIGHT_TYPE")
	if arquivo == "EDGE_WEIGHT_TYPE" {
		arquivo, _ = r.str()
	}

	ewt, _ = r.str()

	x := make([]float64, N+1)
	y := make([]float64, N+1)

	// Alocar matriz 2D
	dist := make([][]float64, N+1)

	for i := 0; i < N+1; i++ {
		dist[i] = make([]float64, N+1)
	}

	if ewt == "EXPLICIT" {

		arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_FORMAT:", "EDGE_WEIGHT_FORMAT")

		var ewf string
		if arquivo == "EDGE_WEIGHT_FORMAT" {
			arquivo, _ = r.str()
		}
		ewf, _ = r.str()

		if ewf == "FUNCTION" {
			fmt.Println("FUNCTION - Nao suportado!")
		} else if ewf == "FULL_MATRIX" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for i := 1; i < N+1; i++ {
				for j := 1; j < N+1; j++ {
					dist[i][j] = r.float()
				}
			}
		} else if ewf == "UPPER_ROW" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for i := 1; i < N; i++ {
				for j := i + 1; j < N+1; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}

			for i := 1; i < N+1; i++ {
				dist[i][i] = 0
			}

		} else if ewf == "LOWER_ROW" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for i := 2; i < N+1; i++ {
				for j := 1; j < i; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}

			for i := 1; i < N+1; i++ {
				dist[i][i] = 0
			}
		} else if ewf == "UPPER_DIAG_ROW" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for i := 1; i < N+1; i++ {
				for j := i; j < N+1; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}
		} else if ewf == "LOWER_DIAG_ROW" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for i := 1; i < N+1; i++ {
				for j := 1; j <= i; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}
		} else if ewf == "UPPER_COL" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for j := 2; j < N+1; j++ {
				for i := 1; i < j; i++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}

			for i := 1; i < N+1; i++ {
				dist[i][i] = 0
			}

		} else if ewf == "LOWER_COL" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for j := 1; j < N; j++ {
				for i := j + 1; i < N+1; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}

			for i := 1; i < N+1; i++ {
				dist[i][i] = 0
			}

		} else if ewf == "UPPER_DIAG_COL" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for j := 1; j < N+1; j++ {
				for i := 1; i <= j; i++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}
		} else if ewf == "LOWER_DIAG_COL" {

			arquivo = r.skipUntil(arquivo, "EDGE_WEIGHT_SECTION")

			// Preencher Matriz Distancia
			for j := 1; j < N+1; j++ {
				for i := j; i < N+1; j++ {
					dist[i][j] = r.float()
					dist[j][i] = dist[i][j]
				}
			}
		}

	} else if ewt == "EUC_2D" {

		arquivo = r.skipUntil(arquivo, "NODE_COORD_SECTION")
		// ler coordenadas
		for i := 1; i < N+1; i++ {
			r.int() // tempCity
			x[i] = r.float()
			y[i] = r.float()
		}

		// Calcular Matriz Distancia (Euclidiana)
		for i := 1; i < N+1; i++ {
			for j := 1; j < N+1; j++ {
				dist[i][j] = math.Floor(CalcDistEuc(x, y, i, j) + 0.5)
			}
		}
	} else if ewt == "EUD_3D" {
		fmt.Println("EUC_3D - Nao suportado!")
	} else if ewt == "MAX_2D" {
		fmt.Println("MAX_2D - Nao suportado!")
	} else if ewt == "MAX_3D" {
		fmt.Println("MAX_3D - Nao Suportado!")
	} else if ewt == "MAN_2D" {
		fmt.Println("MAN_2D - Nao suportado!")
	} else if ewt == "MAN_3D" {
		fmt.Println("MAN_3D - Nao Suportado!")
	} else if ewt == "CEIL_2D" {

		arquivo = r.skipUntil(arquivo, "NODE_COORD_SECTION")
		// ler coordenadas
		for i := 1; i < N+1; i++ {
			r.int() // tempCity
			x[i] = r.float()
			y[i] = r.float()
		}

		// Calcular Matriz Distancia (Euclidiana)
		for i := 1; i < N+1; i++ {
			for j := 1; j < N+1; j++ {
				dist[i][j] = math.Ceil(CalcDistEuc(x, y, i, j))
			}
		}
	} else if ewt == "GEO" {

		arquivo = r.skipUntil(arquivo, "NODE_COORD_SECTION")
		// ler coordenadas
		for i := 1; i < N+1; i++ {
			r.int() // tempCity
			x[i] = r.float()
			y[i] = r.float()
		}

		latitude := make([]float64, N+1)
		longitude := make([]float64, N+1)

		CalcLatLong(x, y, N, latitude, longitude)

		// Calcular Matriz Distancia
		for i := 1; i < N+1; i++ {
			for j := 1; j < N+1; j++ {
				dist[i][j] = CalcDistGeo(latitude, longitude, i, j)
			}
		}

	} else if ewt == "ATT" {

		arquivo = r.skipUntil(arquivo, "NODE_COORD_SECTION")

		// ler coordenadas
		tempX := make([]int, N+1)
		tempY := make([]int, N+1)

		for i := 1; i < N+1; i++ {
			r.int() // tempCity
			tempX[i] = r.int()
			tempY[i] = r.int()
			x[i] = float64(tempX[i])
			y[i] = float64(tempY[i])
		}

		// Calcular Matriz Distancia (Pesudo-Euclidiana)
		for i := 1; i < N+1; i++ {
			for j := 1; j < N+1; j++ {
				dist[i][j] = CalcDistAtt(x, y, i, j)
			}
		}

	} else if ewt == "XRAY1" {
		fmt.Println("XRAY1 - Nao suportado!")
	} else if ewt == "XRAY2" {
		fmt.Println("XRAY2 - Nao suportado!")
	} else if ewt == "SPECIAL" {
		fmt.Println("SPECIAL - Nao suportado!")
	}

	_ = arquivo

	*Dimension = N
	*Mdist = dist
}

func CalcDistEuc(X, Y []float64, I, J int) float64 {
	return math.Sqrt(math.Pow(X[I]-X[J], 2) + math.Pow(Y[I]-Y[J], 2))
}

func CalcDistAtt(X, Y []float64, I, J int) float64 {
	// Calcula Pseudo Distancia Euclidiana
	var rij, tij, dij float64

	rij = math.Sqrt((math.Pow(X[I]-X[J], 2) + math.Pow(Y[I]-Y[J], 2)) / 10)
	tij = math.Floor(rij + 0.5)

	if tij < rij {
		dij = tij + 1
	} else {
		dij = tij
	}

	return dij
}

func CalcLatLong(X, Y []float64, n int, latit, longit []float64) {
	PI := 3.141592
	var min float64
	var deg int

	for i := 1; i < n+1; i++ {
		deg = int(X[i])
		min = X[i] - float64(deg)
		latit[i] = PI * (float64(deg) + 5.0*min/3.0) / 180.0
	}

	for i := 1; i < n+1; i++ {
		deg = int(Y[i])
		min = Y[i] - float64(deg)
		longit[i] = PI * (float64(deg) + 5.0*min/3.0) / 180.0
	}
}

func CalcDistGeo(latit, longit []float64, I, J int) float64 {
	var q1, q2, q3 float64
	RRR := 6378.388

	q1 = math.Cos(longit[I] - longit[J])
	q2 = math.Cos(latit[I] - latit[J])
	q3 = math.Cos(latit[I] + latit[J])

	return float64(int(RRR*math.Acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3)) + 1.0))
}
