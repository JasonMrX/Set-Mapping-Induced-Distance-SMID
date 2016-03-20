clear all;
close all;
files = dir('*.jpg');
Metric = zeros(1, numel(files) / 2);
for i = 1 : 2 : numel(files)
    idx = (i + 1) / 2;
    [mDistPixel, mDistEntropy, mDistDCT, mDistDCTEntropy, success] = ComputeDCTComplexity(files(i).name, files(i + 1).name);
%     title(idx);
    Metric(idx) = mDistDCTEntropy;
    title(mDistDCTEntropy);
%     pause;
end
save Metric Metric;

