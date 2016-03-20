function WriteMetric( SimID, NotSimID, SimMetric, NotSimMetric, SimName, NotSimName, filename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    fout = fopen(filename, 'w');
    fout_csv = fopen('lastPrint.csv', 'w');
    argv = size(SimMetric, 2);
    format = ['%4d %16s %16s %3d   ', repmat('%10.4f ', 1, argv), '\n'];
    format_csv = ['%4d, %16s, %16s, %3d,   ', repmat('%10.4f, ', 1, argv), '\n'];
    for i = 1 : numel(SimID)
        vote = 1;
        fprintf(fout, format, SimID(i), SimName(i).name1, SimName(i).name2, vote, SimMetric(i, :));
        fprintf(fout_csv, format_csv, SimID(i), SimName(i).name1, SimName(i).name2, vote, SimMetric(i, :));
    end
    for i = 1 : numel(NotSimID)
        vote = 2;
        fprintf(fout, format, NotSimID(i), NotSimName(i).name1, NotSimName(i).name2, vote, NotSimMetric(i, :));
        fprintf(fout_csv, format_csv, NotSimID(i), NotSimName(i).name1, NotSimName(i).name2, vote, NotSimMetric(i, :));
    end
    fclose(fout);
    fclose(fout_csv);
end

