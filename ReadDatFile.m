function [ImagePair, Metric] = ReadDatFile(filename, argc)
    fin = fopen(filename, 'r');
    len = 1 ;
    pos1 = 1 ;
    while (len > 0)
        [imageid, len] = fscanf(fin, '%d ', 1);
        if len > 0
            name1 = fscanf(fin, '%s ', 1);
            name2 = fscanf(fin, '%s ', 1);
            type = fscanf(fin, '%d ', 1);
            ImagePair(pos1).imageid = imageid;

            ImagePair(pos1).name1 = name1;
            ImagePair(pos1).name2 = name2;
            ImagePair(pos1).vote = type;
            Metric(pos1, 1 : argc) = fscanf(fin, '%f ', argc);
            pos1 = pos1 + 1 ;
        end
    end
    fclose(fin);
end