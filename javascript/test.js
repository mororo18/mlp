

N = 3000

m_a = []
for (var i = 0; i < N; i++) {
    m_a[i] = Array(N);
    for (var j = 0; j < N; j++) {
        m_a[i][j] = Array(3);
        m_a[i][j][0] = 0.0;
        m_a[i][j][1] = 0.0;
        m_a[i][j][2] = 0.0;

    }
}

m_b = []
for (var i = 0; i < N; i++) {
    m_b[i] = Array(N*3)
    for (var j = 0; j < N; j++) {
        m_b[i*N*3 + j*3 + 0] = 0.0;
        m_b[i*N*3 + j*3 + 1] = 0.0;
        m_b[i*N*3 + j*3 + 2] = 0.0;
    }
}

index = []
for (var i = 0; i < N; i++) {
    index[i] = parseInt(Math.random() * N)
}

var start = new Date();
for (var i = 0; i < N; i++) {
    for (var j = 0; j < N; j++) {
        var x = m_a[index[i]][index[j]][0];
        var y = m_a[index[i]][index[j]][1];
        var z = m_a[index[i]][index[j]][2];
    }
}

var end = new Date();
m_a.length = 0;
console.log("TIME  A: ", (end-start)/1000);

var to_1D = function (x, y, z, n) {
    return (x*n + y)*3 + z;
}



start = new Date();
for (var i = 0; i < N; i++) {
    for (var j = 0; j < N; j++) {
      //var x = m_b[(index[i]*N + index[j])*3 + 0];
      //var y = m_b[(index[i]*N + index[j])*3 + 1];
      //var z = m_b[(index[i]*N + index[j])*3 + 2];
        var x = m_b[to_1D(index[i], index[j],  0, N)];
        var y = m_b[to_1D(index[i], index[j],  1, N)];
        var z = m_b[to_1D(index[i], index[j],  2, N)];
    }
}
end = new Date();
console.log("TIME  B: ", (end-start)/1000);

function to_1D(i, j, k, n) {

}
