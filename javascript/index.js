let hello = "HELLO WORLD!!";

/*
var getFileBlob = function (url, cb) {
    var xhr = new XMLHttpRequest();
    xhr.open("GET", url);
    xhr.responseType = "blob";
    xhr.addEventListener('load', function() {
        cb(xhr.response);
    });
    xhr.send();
};

var blobToFile = function (blob, name) {
    blob.lastModifiedDate = new Date();
    blob.name = name;
    return blob;
};

var getFileObject = function(filePathOrUrl, cb) {
    getFileBlob(filePathOrUrl, function (blob) {
        cb(blobToFile(blob, 'distence_matrix'));
    });
};

getFileObject('file:///home/victor//repos/mlp/distance_matrix', function (fileObject) {
    console.log(fileObject);
}); 
*/

function printFile(file) {
    const reader = new FileReader();
    reader.onload = function(evt) {
        console.log(evt.target.result);
    };
    reader.readAsText(file);
}

//printFile('file:///home/victor//repos/mlp/distance_matrix');

function getText(){
    /*
    var request = new XMLHttpRequest();
    request.open('GET', 'file:///home/victor//repos/mlp/distance_matrix');
    request.send(null);
    request.onreadystatechange = function () {
        if (request.readyState === 4 && request.status === 200) {
            var type = request.getResponseHeader('Content-Type');
            if (type.indexOf("text") !== 1) {
                return request.responseText;
            }
        }
    }
*/    
    document.getElementById('inputfile')
        .addEventListener('change', function() {
            var file;

            var fr=new FileReader();
            fr.onload=function(){
                document.getElementById('output')
                    .textContent=fr.result;
                file = fr.result;
            }

            fr.readAsBinaryString(this.files[0]);
            //console.log(file.length);
        })

    var fs = require("fs");

}
console.log(getText());
