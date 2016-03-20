function [ rocCurve ] = ComputeROC( SimMetric, NotSimMetric, type, resolution )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 2
        type = 1;
    end
    if type == 2
        tmp = SimMetric;
        SimMetric = NotSimMetric;
        NotSimMetric = tmp;
    end
    low = 0;
    high = max([SimMetric; NotSimMetric]);
    if nargin <= 3
        resolution = (high - low) / 1000;   
    end
    range = low: resolution: high;
    rocCurve = [zeros(size(repmat(range, 2, 1))), [0; 0]];
    for i = 1 : numel(range)
        rocCurve(1, i + 1) = sum(NotSimMetric < range(i)) / numel(NotSimMetric);
        rocCurve(2, i + 1) = sum(SimMetric < range(i)) / numel(SimMetric);
    end
    if type == 2
        rocCurve = 1 - rocCurve;
        rocCurve = rocCurve([2 1], :);
    end
end

