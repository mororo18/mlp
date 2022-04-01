var fs = require("fs");
var text = fs.readFileSync("../distance_matrix").toString('utf-8');
var file = text.split("\n");

function matrix_fill(matrix, n) {
    for (var i = 0; i < n; i++) {
        matrix[i] = new Array(n);
    }
}

module.exports = {
    info_load: function(c, rnd) {
        let index = 0;
        console.log(file[index]);

        let dimension = parseInt(file[index++]);
        matrix_fill(c, dimension);

        for (var i = 0; i < dimension; i++) {
            c[i][i] = 0.0;
            cost = file[index].split(" ");
            for (var j = i+1; j < dimension; j++) {
                c[i][j] = parseFloat(cost[j-(i+1)]);
                c[j][i] = parseFloat(cost[j-(i+1)]);
            }

            index++;
        }

        index++;
        index++;
        console.log(file[index]);
        let rnd_size = parseInt(file[index++]);
        rnd = new Array(rnd_size);

        for (var i = 0; i < rnd_size; i++) {
            rnd[i] = file[index++];
        }

        console.log(rnd);

        /*
        for (var i = -1; i+1 < file.length; i++) {
            //console.log(file[i+1]);
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
        */
        return dimension;
    }

}
