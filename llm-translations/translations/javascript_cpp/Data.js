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
        let index = 0;

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
        let rnd_size = parseInt(file[index++]);
        var rnd = new Array(rnd_size);

        for (var i = 0; i < rnd_size; i++) {
            rnd[i] = parseInt(file[index++]);
        }

        return {
            dimension: dimension, 
            rnd : rnd
        };
    }

}
