function [dim, cost] = Data
    txt = fileread("../distance_matrix");
    lines = strsplit(txt, "\n");
    typeinfo(lines);
    c = [];
    dimension = NaN;
    i = 0;
    for line = lines
        str = line{1}(:);
        typeinfo(str);
        if (!isempty(strfind(str, "EOF")))
            line;
            break
        end

        if (i == 0)
            v = mat2str(str2double(str));
            v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
            dimension = v;
            c = zeros(dimension, dimension);
            i = i + 1;
            continue;
        end

        j = i + 1;
        while (true)

            pos = strfind(str, " ");
            if (ismatrix(pos) && isempty(pos))
                %"yes"
                break;
            end
            pos = pos(1, 1);
            v = mat2str(str2double(str(1:pos-1)));
            v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
            c(i, j) = v;
            typeinfo(v);
            str = str(pos+1:end);
            %break;
            j = j + 1;
        end

        i = i + 1;
    end

    dim = dimension;
    cost = c;
end
