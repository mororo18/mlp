import sys
import math


def calc_dist_euc(x, y, i, j):
    return math.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2)


def calc_dist_att(x, y, i, j):
    rij = math.sqrt(((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2) / 10)
    tij = math.floor(rij + 0.5)

    if tij < rij:
        dij = tij + 1
    else:
        dij = tij

    return dij


def calc_lat_long(x, y, n):
    pi = 3.141592
    latit = [0.0] * (n + 1)
    longit = [0.0] * (n + 1)

    for i in range(1, n + 1):
        deg = int(x[i])
        minutes = x[i] - deg
        latit[i] = pi * (deg + 5.0 * minutes / 3.0) / 180.0

    for i in range(1, n + 1):
        deg = int(y[i])
        minutes = y[i] - deg
        longit[i] = pi * (deg + 5.0 * minutes / 3.0) / 180.0

    return latit, longit


def calc_dist_geo(latit, longit, i, j):
    rrr = 6378.388

    q1 = math.cos(longit[i] - longit[j])
    q2 = math.cos(latit[i] - latit[j])
    q3 = math.cos(latit[i] + latit[j])

    return int(rrr * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)


def _read_tokens(f):
    for line in f:
        for tok in line.split():
            yield tok


def read_data(argv):
    if len(argv) < 2:
        print("\nFaltando parametros\n")
        print(" ./exec [Instancia] ")
        sys.exit(1)

    if len(argv) > 2:
        print("\nMuitos parametros\n")
        print(" ./exec [Instancia] ")
        sys.exit(1)

    instancia = argv[1]

    try:
        f = open(instancia, "r")
    except IOError:
        print("arquivo nao pode ser aberto\n")
        sys.exit(1)

    tokens = _read_tokens(f)

    def next_tok():
        return next(tokens)

    arquivo = ""
    while arquivo != "DIMENSION:" and arquivo != "DIMENSION":
        arquivo = next_tok()

    if arquivo == "DIMENSION":
        arquivo = next_tok()

    n = int(next_tok())

    while arquivo != "EDGE_WEIGHT_TYPE:" and arquivo != "EDGE_WEIGHT_TYPE":
        arquivo = next_tok()
    if arquivo == "EDGE_WEIGHT_TYPE":
        arquivo = next_tok()

    ewt = next_tok()

    x = [0.0] * (n + 1)
    y = [0.0] * (n + 1)

    dist = [[0.0] * (n + 1) for _ in range(n + 1)]

    if ewt == "EXPLICIT":
        while arquivo != "EDGE_WEIGHT_FORMAT:" and arquivo != "EDGE_WEIGHT_FORMAT":
            arquivo = next_tok()

        if arquivo == "EDGE_WEIGHT_FORMAT":
            arquivo = next_tok()
        ewf = next_tok()

        if ewf == "FUNCTION":
            print("FUNCTION - Nao suportado!")

        elif ewf == "FULL_MATRIX":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for i in range(1, n + 1):
                for j in range(1, n + 1):
                    dist[i][j] = float(next_tok())

        elif ewf == "UPPER_ROW":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for i in range(1, n):
                for j in range(i + 1, n + 1):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

            for i in range(1, n + 1):
                dist[i][i] = 0

        elif ewf == "LOWER_ROW":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for i in range(2, n + 1):
                for j in range(1, i):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

            for i in range(1, n + 1):
                dist[i][i] = 0

        elif ewf == "UPPER_DIAG_ROW":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for i in range(1, n + 1):
                for j in range(i, n + 1):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

        elif ewf == "LOWER_DIAG_ROW":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for i in range(1, n + 1):
                for j in range(1, i + 1):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

        elif ewf == "UPPER_COL":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for j in range(2, n + 1):
                for i in range(1, j):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

            for i in range(1, n + 1):
                dist[i][i] = 0

        elif ewf == "LOWER_COL":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for j in range(1, n):
                i = j + 1
                while i < n + 1:
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]
                    j += 1

            for i in range(1, n + 1):
                dist[i][i] = 0

        elif ewf == "UPPER_DIAG_COL":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for j in range(1, n + 1):
                for i in range(1, j + 1):
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]

        elif ewf == "LOWER_DIAG_COL":
            while arquivo != "EDGE_WEIGHT_SECTION":
                arquivo = next_tok()

            for j in range(1, n + 1):
                i = j
                while i < n + 1:
                    dist[i][j] = float(next_tok())
                    dist[j][i] = dist[i][j]
                    j += 1

    elif ewt == "EUC_2D":
        while arquivo != "NODE_COORD_SECTION":
            arquivo = next_tok()
        for i in range(1, n + 1):
            next_tok()
            x[i] = float(next_tok())
            y[i] = float(next_tok())

        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dist[i][j] = math.floor(calc_dist_euc(x, y, i, j) + 0.5)

    elif ewt == "EUD_3D":
        print("EUC_3D - Nao suportado!")

    elif ewt == "MAX_2D":
        print("MAX_2D - Nao suportado!")

    elif ewt == "MAX_3D":
        print("MAX_3D - Nao Suportado!")

    elif ewt == "MAN_2D":
        print("MAN_2D - Nao suportado!")

    elif ewt == "MAN_3D":
        print("MAN_3D - Nao Suportado!")

    elif ewt == "CEIL_2D":
        while arquivo != "NODE_COORD_SECTION":
            arquivo = next_tok()
        for i in range(1, n + 1):
            next_tok()
            x[i] = float(next_tok())
            y[i] = float(next_tok())

        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dist[i][j] = math.ceil(calc_dist_euc(x, y, i, j))

    elif ewt == "GEO":
        while arquivo != "NODE_COORD_SECTION":
            arquivo = next_tok()
        for i in range(1, n + 1):
            next_tok()
            x[i] = float(next_tok())
            y[i] = float(next_tok())

        latitude, longitude = calc_lat_long(x, y, n)

        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dist[i][j] = calc_dist_geo(latitude, longitude, i, j)

    elif ewt == "ATT":
        while arquivo != "NODE_COORD_SECTION":
            arquivo = next_tok()

        temp_x = [0] * (n + 1)
        temp_y = [0] * (n + 1)
        for i in range(1, n + 1):
            next_tok()
            temp_x[i] = int(next_tok())
            temp_y[i] = int(next_tok())
            x[i] = temp_x[i]
            y[i] = temp_y[i]

        for i in range(1, n + 1):
            for j in range(1, n + 1):
                dist[i][j] = calc_dist_att(x, y, i, j)

    elif ewt == "XRAY1":
        print("XRAY1 - Nao suportado!")

    elif ewt == "XRAY2":
        print("XRAY2 - Nao suportado!")

    elif ewt == "SPECIAL":
        print("SPECIAL - Nao suportado!")

    f.close()

    return n, dist
