function [dim, cost, ret_rnd] = Data
    txt = fileread("../distance_matrix");
    lines = strsplit(txt, "\n");
    c = [];
    dimension = NaN;
    i = 0;

    index = 1;
    dimension = str2double(lines{index});
    index = index + 1;

    c = zeros(dimension, dimension);

    for i = 1:dimension
        j = i + 1;

        str = lines(index);
        str = str{1};
        class(str);
        for j = i+1:dimension
            class(str);
            pos = strfind(str, " ");
            if (ismatrix(pos) && isempty(pos))
                break;
            end
            pos = pos(1, 1);
            v = mat2str(str2double(str(1:pos-1)));
            
            a = size(v);
            v = str2double(v);
            c(i, j) = v;
            c(j, i) = v;
            class(v);
            str = str(pos+1:end);
        end
        index = index + 1;
    end
    
    index = index + 2;
    rnd_size = str2double(lines{index});
    rnd = [];
    for i = 1:rnd_size
        v = str2double(lines{index + i});
        rnd(i) = v;
    end
    
    dim = dimension;
    cost = c;
    ret_rnd = rnd;
end
