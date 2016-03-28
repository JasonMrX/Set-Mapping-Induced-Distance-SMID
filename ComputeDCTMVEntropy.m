function [ EntropyArray ] = ComputeDCTMVEntropy( filename, width, height )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    a = dlmread(filename);
    for i = 1 : 6
        py = reshape(a(:, (i - 1) * 2 + 1), height, width)';
        px = reshape(a(:, i * 2), height, width)';
        ey = py - spatialPredict(py);
        ex = px - spatialPredict(px);
        EntropyArray((i - 1) * 4 + 1 : i * 4) = [calEntropy(py(:)), calEntropy(ey(:)), calEntropy(px(:)), calEntropy(ex(:))];
    end

end

