clear all;
close all;
files = dir('*.jpg');
load Metric;
Entropy = Metric;
load Metric_DCTDist;
DCTDist = Metric;
for i = 1 : 2 : numel(files)
    idx = (i + 1) / 2;
    [mDistPixel, mDistEntropy, mDistDCT, mDistDCTEntropy, success] = ComputeDCTComplexity(files(i).name, files(i + 1).name);
    pause;
    title([idx, DCTDist(idx)]);
    pause; 
end

