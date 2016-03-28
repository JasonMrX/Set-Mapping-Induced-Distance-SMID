function [ pred ] = spatialPredict( org )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    pred = org;
    pred(1, 2 : end) = org(1, 2 : end) - org(1, 1 : end - 1);
    pred(2 : end, 1) = org(2 : end, 1) - org(1 : end - 1, 1);
    A = org(1 : end - 1, 2 : end);
    B = org(2 : end, 1 : end - 1);
    C = org(2 : end, 2 : end);
    D = pred(2 : end, 2 : end);
    
    minAB = min(A, B);
    maxAB = max(A, B);
    AB_C = A + B - C;
    
    idxMat = (C >= maxAB);
    D(idxMat) = minAB(idxMat);
    
    idxMat = (C <= minAB);
    D(idxMat) = maxAB(idxMat);
    
    idxMat = ((C < maxAB) .* (C > minAB)) == 1;
    D(idxMat) = AB_C(idxMat);
    
    pred(2 : end, 2 : end) = D;

end

