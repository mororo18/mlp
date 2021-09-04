var fs = require("fs");
var text = fs.readFileSync("../distance_matrix").toString('utf-8');
var file = text.split("\n");

function matrix_fill(matrix, n) {
    for (var i = 0; i < n; i++) {
        matrix[i] = new Array(n);
    }
}

module.exports = {
    info_load: function(c) {
        for (var i = -1; i+1 < file.length; i++) {
            console.log(file[i+1]);
            let l = i+1;
            if (i == -1) {
                var dimension = parseInt(file[l]);
                matrix_fill(c, dimension);
                continue;
            }

            if (file[l] == "EOF") {
                break;
            }

            let j = i + 1;
            cost = file[l].split(" ");
            c[i][i] = 0.0;
            for (k in cost) {
                if (j >= dimension)
                    break;
                c[i][j] = parseFloat(cost[k]);
                c[j][i] = parseFloat(cost[k]);
                j++;
            }
        }
        c[dimension-1][dimension-1] = 0.0;
        return dimension;
    }

}
