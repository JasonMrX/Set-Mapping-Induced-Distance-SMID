function [entropy] = calEntropy(data)
    h = hist(data, min(data):max(data));
    h = h(h ~= 0);
    h = h / sum(h);
    entropy = sum(-h .* log2(h));
end