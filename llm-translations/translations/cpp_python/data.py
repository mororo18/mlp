def matrix_allocate(sz):
    return [[0.0] * sz for _ in range(sz)]


def load_data():
    rnd = []
    dimension = 0
    matrix = None

    with open("../distance_matrix", "r") as file:
        lines = file.read().splitlines()

    if not lines:
        return dimension, matrix, rnd

    dimension = int(lines[0])
    matrix = matrix_allocate(dimension)

    flag = False
    rnd_size = -1
    i = 0
    for line in lines[1:]:
        j = i + 1
        while " " in line:
            pos = line.find(" ")
            value = float(line[:pos])
            line = line[pos + 1:]
            matrix[i][j] = value
            matrix[j][i] = value
            j += 1

        i += 1

        if flag:
            if rnd_size < 0:
                rnd_size = int(line)
            else:
                rnd.append(int(line))

        if line == "RND":
            flag = True

    return dimension, matrix, rnd


class Data:
    def __init__(self):
        self.dimen = 0
        self.cost = None
        self.rnd = []

    def load(self):
        self.dimen, self.cost, self.rnd = load_data()

    def get_rnd_vec(self):
        return self.rnd

    def get_dimen(self):
        return self.dimen

    def get_cost_ptr(self):
        return self.cost
