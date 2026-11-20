import Foundation

func calcDistEuc(_ X: [Double], _ Y: [Double], _ I: Int, _ J: Int) -> Double {
    return (pow(X[I] - X[J], 2) + pow(Y[I] - Y[J], 2)).squareRoot()
}

func calcDistAtt(_ X: [Double], _ Y: [Double], _ I: Int, _ J: Int) -> Double {
    let rij = ((pow(X[I] - X[J], 2) + pow(Y[I] - Y[J], 2)) / 10).squareRoot()
    let tij = (rij + 0.5).rounded(.down)

    let dij: Double
    if tij < rij {
        dij = tij + 1
    } else {
        dij = tij
    }

    return dij
}

func calcLatLong(_ X: [Double], _ Y: [Double], _ n: Int, _ latit: inout [Double], _ longit: inout [Double]) {
    let PI = 3.141592

    for i in 1...n {
        let deg = Double(Int(X[i]))
        let min = X[i] - deg
        latit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0
    }

    for i in 1...n {
        let deg = Double(Int(Y[i]))
        let min = Y[i] - deg
        longit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0
    }
}

func calcDistGeo(_ latit: [Double], _ longit: [Double], _ I: Int, _ J: Int) -> Double {
    let RRR = 6378.388

    let q1 = cos(longit[I] - longit[J])
    let q2 = cos(latit[I] - latit[J])
    let q3 = cos(latit[I] + latit[J])

    return Double(Int(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0))
}

private struct TokenReader {
    private let tokens: [String]
    private var index = 0

    init(_ text: String) {
        tokens = text.split(whereSeparator: { $0 == " " || $0 == "\n" || $0 == "\t" || $0 == "\r" }).map(String.init)
    }

    mutating func next() -> String {
        guard index < tokens.count else { return "" }
        let t = tokens[index]
        index += 1
        return t
    }

    mutating func nextInt() -> Int {
        return Int(next()) ?? 0
    }

    mutating func nextDouble() -> Double {
        return Double(next()) ?? 0.0
    }
}

func readData(_ argc: Int, _ argv: [String], _ Dimension: inout Int, _ Mdist: inout [[Double]]) {

    if argc < 2 {
        print("\nFaltando parametros")
        print(" ./exec [Instancia] ")
        exit(1)
    }

    if argc > 2 {
        print("\nMuitos parametros")
        print(" ./exec [Instancia] ")
        exit(1)
    }

    var N = 0
    var arquivo = ""
    var ewt = ""

    let instancia = argv[1]

    guard let content = try? String(contentsOfFile: instancia, encoding: .utf8) else {
        print("arquivo nao pode ser aberto")
        exit(1)
    }

    var reader = TokenReader(content)

    while arquivo != "DIMENSION:" && arquivo != "DIMENSION" {
        arquivo = reader.next()
    }

    if arquivo == "DIMENSION" { arquivo = reader.next() }

    N = reader.nextInt()

    while arquivo != "EDGE_WEIGHT_TYPE:" && arquivo != "EDGE_WEIGHT_TYPE" {
        arquivo = reader.next()
    }
    if arquivo == "EDGE_WEIGHT_TYPE" { arquivo = reader.next() }

    ewt = reader.next()

    var x = Array(repeating: 0.0, count: N + 1)
    var y = Array(repeating: 0.0, count: N + 1)

    var dist = Array(repeating: Array(repeating: 0.0, count: N + 1), count: N + 1)

    if ewt == "EXPLICIT" {

        while arquivo != "EDGE_WEIGHT_FORMAT:" && arquivo != "EDGE_WEIGHT_FORMAT" {
            arquivo = reader.next()
        }

        var ewf = ""
        if arquivo == "EDGE_WEIGHT_FORMAT" { arquivo = reader.next() }
        ewf = reader.next()

        if ewf == "FUNCTION" {
            print("FUNCTION - Nao suportado!")
        }

        else if ewf == "FULL_MATRIX" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for i in 1..<(N + 1) {
                for j in 1..<(N + 1) {
                    dist[i][j] = reader.nextDouble()
                }
            }
        }

        else if ewf == "UPPER_ROW" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for i in 1..<N {
                for j in (i + 1)..<(N + 1) {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }

            for i in 1..<(N + 1) {
                dist[i][i] = 0
            }
        }

        else if ewf == "LOWER_ROW" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for i in 2..<(N + 1) {
                for j in 1..<i {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }

            for i in 1..<(N + 1) {
                dist[i][i] = 0
            }
        }

        else if ewf == "UPPER_DIAG_ROW" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for i in 1..<(N + 1) {
                for j in i..<(N + 1) {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }
        }

        else if ewf == "LOWER_DIAG_ROW" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for i in 1..<(N + 1) {
                for j in 1...i {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }
        }

        else if ewf == "UPPER_COL" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for j in 2..<(N + 1) {
                for i in 1..<j {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }

            for i in 1..<(N + 1) {
                dist[i][i] = 0
            }
        }

        else if ewf == "LOWER_COL" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            // NOTE: mirrors the original C++'s non-incrementing outer index bug (j never advances)
            for j in 1..<N {
                var i = j + 1
                while i < N + 1 {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                    i = j + 1
                }
            }

            for i in 1..<(N + 1) {
                dist[i][i] = 0
            }
        }

        else if ewf == "UPPER_DIAG_COL" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            for j in 1..<(N + 1) {
                for i in 1...j {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                }
            }
        }

        else if ewf == "LOWER_DIAG_COL" {

            while arquivo != "EDGE_WEIGHT_SECTION" {
                arquivo = reader.next()
            }

            // NOTE: mirrors the original C++'s non-incrementing outer index bug (j never advances)
            for j in 1..<(N + 1) {
                var i = j
                while i < N + 1 {
                    dist[i][j] = reader.nextDouble()
                    dist[j][i] = dist[i][j]
                    i = j
                }
            }
        }

    }

    else if ewt == "EUC_2D" {

        while arquivo != "NODE_COORD_SECTION" {
            arquivo = reader.next()
        }

        for i in 1..<(N + 1) {
            _ = reader.nextInt()
            x[i] = reader.nextDouble()
            y[i] = reader.nextDouble()
        }

        for i in 1..<(N + 1) {
            for j in 1..<(N + 1) {
                dist[i][j] = (calcDistEuc(x, y, i, j) + 0.5).rounded(.down)
            }
        }
    }

    else if ewt == "EUD_3D" {
        print("EUC_3D - Nao suportado!")
    }

    else if ewt == "MAX_2D" {
        print("MAX_2D - Nao suportado!")
    }

    else if ewt == "MAX_3D" {
        print("MAX_3D - Nao Suportado!")
    }

    else if ewt == "MAN_2D" {
        print("MAN_2D - Nao suportado!")
    }

    else if ewt == "MAN_3D" {
        print("MAN_3D - Nao Suportado!")
    }

    else if ewt == "CEIL_2D" {

        while arquivo != "NODE_COORD_SECTION" {
            arquivo = reader.next()
        }

        for i in 1..<(N + 1) {
            _ = reader.nextInt()
            x[i] = reader.nextDouble()
            y[i] = reader.nextDouble()
        }

        for i in 1..<(N + 1) {
            for j in 1..<(N + 1) {
                dist[i][j] = calcDistEuc(x, y, i, j).rounded(.up)
            }
        }
    }

    else if ewt == "GEO" {

        while arquivo != "NODE_COORD_SECTION" {
            arquivo = reader.next()
        }

        for i in 1..<(N + 1) {
            _ = reader.nextInt()
            x[i] = reader.nextDouble()
            y[i] = reader.nextDouble()
        }

        var latitude = Array(repeating: 0.0, count: N + 1)
        var longitude = Array(repeating: 0.0, count: N + 1)

        calcLatLong(x, y, N, &latitude, &longitude)

        for i in 1..<(N + 1) {
            for j in 1..<(N + 1) {
                dist[i][j] = calcDistGeo(latitude, longitude, i, j)
            }
        }
    }

    else if ewt == "ATT" {

        while arquivo != "NODE_COORD_SECTION" {
            arquivo = reader.next()
        }

        var tempX = Array(repeating: 0, count: N + 1)
        var tempY = Array(repeating: 0, count: N + 1)

        for i in 1..<(N + 1) {
            _ = reader.nextInt()
            tempX[i] = reader.nextInt()
            tempY[i] = reader.nextInt()
            x[i] = Double(tempX[i])
            y[i] = Double(tempY[i])
        }

        for i in 1..<(N + 1) {
            for j in 1..<(N + 1) {
                dist[i][j] = calcDistAtt(x, y, i, j)
            }
        }
    }

    else if ewt == "XRAY1" {
        print("XRAY1 - Nao suportado!")
    }

    else if ewt == "XRAY2" {
        print("XRAY2 - Nao suportado!")
    }

    else if ewt == "SPECIAL" {
        print("SPECIAL - Nao suportado!")
    }

    Dimension = N
    Mdist = dist
}
