'use strict';

const fs = require('fs');

function matrix_allocate(sz) {
    const ptr = new Array(sz);

    for (let i = 0; i < sz; i++) {
        ptr[i] = new Float64Array(sz);
    }
    return ptr;
}

// Reads "../distance_matrix" the same way the C++ loadData() does:
// first line is the dimension, then the upper triangle of the cost matrix
// (one row per line, values separated by spaces), then a "RND" marker
// followed by the amount of random numbers and the numbers themselves.
function loadData(rnd) {
    let dimension;
    let matrix = null;

    let content;
    try {
        content = fs.readFileSync('../distance_matrix', 'utf8');
    } catch (e) {
        return { dimension: undefined, matrix: null };
    }

    const lines = content.split('\n');
    let lineIdx = 0;

    let line = lines[lineIdx++];
    dimension = parseInt(line, 10);
    matrix = matrix_allocate(dimension);

    let flag = false;
    let rnd_size = -1;
    let i = 0;
    while (lineIdx < lines.length) {
        line = lines[lineIdx++];
        if (line.length > 0 && line.charCodeAt(line.length - 1) === 13) {
            line = line.slice(0, -1);
        }
        // getline() stops before EOF; an empty trailing chunk is not a line
        if (line === '' && lineIdx >= lines.length) break;

        let j = i + 1;
        let pos;
        while ((pos = line.indexOf(' ')) !== -1) {
            const value = parseFloat(line.substring(0, pos));
            line = line.substring(pos + 1);
            matrix[i][j] = value;
            matrix[j][i] = value;
            j++;
        }

        i++;

        if (flag) {

            if (rnd_size < 0) {
                rnd_size = parseInt(line, 10);
            } else {
                rnd.push(parseInt(line, 10));
            }
        }

        if (line === 'RND') flag = true;
    }

    return { dimension, matrix };
}

module.exports = { loadData };
