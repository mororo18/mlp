function [dim, cost, ret_rnd] = Data
    txt = fileread("../distance_matrix");
    lines = strsplit(txt, "\n");
    %typeinfo(lines);
    c = [];
    dimension = NaN;
    i = 0;

    index = 1;
    dimension = mat2str(str2double(lines{index}(:)));
    dimension = str2double(dimension);
    %dimension = str2double(strrep(substr(dimension, 2, numel(dimension)-2), ';', ''));

    index = index + 1;

    c = zeros(dimension, dimension);

    for i = 1:dimension
        j = i + 1;

        str = lines(index);
        str = str{1};
        class(str);
        %str = str(:);
        for j = i+1:dimension
            %lines(index){1};
            %typeinfo(str);
            class(str);
            pos = strfind(str, " ");
            if (ismatrix(pos) && isempty(pos))
                %"yes"
                break;
            end
            pos = pos(1, 1);
            v = mat2str(str2double(str(1:pos-1)));
            
            a = size(v);
            %a(1,2)
            %if (a(1,2) > 1)
            %    v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
            %else
                v = str2double(v);
            %end
            c(i, j) = v;
            c(j, i) = v;
            class(v);
            str = str(pos+1:end);
        end
        index = index + 1;
    end
    
    index = index + 2;
    rnd_size = str2double(mat2str(str2double(lines{index}(:))));
    %rnd_size = str2double(strrep(substr(rnd_size, 2, numel(rnd_size)-2), ';', ''));

    rnd = [];
    for i = 1:rnd_size
        v = mat2str(str2double(lines{index + i}(:)));
        lines{index + i};
        class(v);
        %if (numel(v) > 1)
        %    v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
        %else 
            v = str2double(v);
        %end
        rnd(i) = v;
    end
    

   %typeinfo(lines)
   %typeinfo(lines{2})
   %for line = lines
   %    str = line{1}(:);
   %    typeinfo(str);
   %    if (!isempty(strfind(str, "EOF")))
   %        line;
   %        break
   %    end

   %    if (i == 0)
   %        v = mat2str(str2double(str));
   %        v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
   %        dimension = v;
   %        c = zeros(dimension, dimension);
   %        i = i + 1;
   %        continue;
   %    end

   %    j = i + 1;
   %    while (true)

   %        pos = strfind(str, " ");
   %        if (ismatrix(pos) && isempty(pos))
   %            %"yes"
   %            break;
   %        end
   %        pos = pos(1, 1);
   %        v = mat2str(str2double(str(1:pos-1)));

   %        if (size(v)(1,2) > 1)
   %            v = str2double(strrep(substr(v, 2, numel(v)-2), ';', ''));
   %        else
   %            v = str2double(v);
   %        end
   %        c(i, j) = v;
   %        c(j, i) = v;
   %        typeinfo(v);
   %        str = str(pos+1:end);
   %        %break;
   %        j = j + 1;
   %    end

   %    i = i + 1;
   %end

    dim = dimension;
    cost = c;
    ret_rnd = rnd;
end
