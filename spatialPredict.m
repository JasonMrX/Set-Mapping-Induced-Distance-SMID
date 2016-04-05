function [ pred, e ] = spatialPredict( org, mask )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if nargin > 2
        error('too many inputs')
    elseif nargin == 1
        mask = ones(size(org));
    end
    
    pred = zeros(size(org));
    pred(1, 2 : end) = org(1, 2 : end) - org(1, 1 : end - 1);
    pred(2 : end, 1) = org(2 : end, 1) - org(1 : end - 1, 1);
    A = org(1 : end - 1, 2 : end);
    AM = mask(1 : end - 1, 2 : end);
    B = org(2 : end, 1 : end - 1);
    BM = mask(2 : end, 1 : end - 1);
    C = org(2 : end, 2 : end);
    CM = mask(2 : end, 2 : end);
    D = pred(2 : end, 2 : end);
    
    if nargin == 1
        minAB = min(A, B);
        maxAB = max(A, B);
        AB_C = A + B - C;

        idxMat = (C >= maxAB);
        D(idxMat) = minAB(idxMat);

        idxMat = (C <= minAB);
        D(idxMat) = maxAB(idxMat);

        idxMat = ((C < maxAB) .* (C > minAB)) == 1;
        D(idxMat) = AB_C(idxMat);
    else
        D(CM == 1) = C(CM == 1);
        D(BM == 1) = B(BM == 1);
        D(AM == 1) = A(AM == 1);
    end
    
    pred(2 : end, 2 : end) = D;
    e = pred - org;
end

