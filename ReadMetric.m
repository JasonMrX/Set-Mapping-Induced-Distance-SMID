function [SimID, NotSimID, SimMetric, NotSimMetric, SimName, NotSimName] = ReadMetric(filename, argv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    numData = 2668;
    fin = fopen(filename, 'r');
    len = 1 ;
    pos1 = 1 ;
    pos2 = 1 ;
    SimID = zeros(numData, 1);
    NotSimID = zeros(numData, 1);
    SimMetric = zeros(numData, argv);
    NotSimMetric = zeros(numData, argv);
    SimName(numData).name1 = '';
    SimName(numData).name2 = '';
    NotSimName(numData).name1 = '';
    NotSimName(numData).name2 = '';
    while(len > 0)
        [imageid, len]= fscanf(fin, '%d ', 1);
        if len > 0
            name1 = fscanf(fin, '%s ', 1);
            name2 = fscanf(fin, '%s ', 1);
            type = fscanf(fin, '%d ', 1);
            if type == 1
                SimMetric(pos1, 1 : argv) = fscanf(fin, '%f ', argv);
                SimID(pos1) = imageid;
                SimName(pos1).name1 = name1 ;
                SimName(pos1).name2 = name2 ;
                pos1 = pos1 + 1 ;
            end
            if type == 2
                NotSimMetric(pos2, 1 : argv) = fscanf(fin, '%f ', argv);
                NotSimID(pos2) = imageid;
                NotSimName(pos2).name1 = name1 ;
                NotSimName(pos2).name2 = name2 ;
                pos2 = pos2 + 1 ;
            end
        end
    end
    SimID = SimID(1 : pos1 - 1);
    SimMetric = SimMetric(1 : pos1 - 1, :);
    SimName = SimName(1 : pos1 - 1);
    NotSimID = NotSimID(1 : pos2 - 1);
    NotSimMetric = NotSimMetric(1 : pos2 - 1, :);
    NotSimName = NotSimName(1 : pos2 - 1);
    fclose(fin);
end

