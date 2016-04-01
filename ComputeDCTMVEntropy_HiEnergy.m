function ComputeDCTMVEntropy_HiEnergy( mvfilename, filename1, filename2 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    sizef = 264*264;
   [ X1oriL] = ReadResizeGrayImage(filename1, sizef) ;
   [ X2oriL] = ReadResizeGrayImage(filename2, sizef) ;
    
%    close all;
%    figure
%    imshow(uint8(X1oriL));
%    figure
%    imshow(uint8(X2oriL));
   
   % watch the same size
   [rowL, colL] = size(X1oriL);
   [row2, col2] = size(X2oriL);
   if rowL ~= row2 || colL ~= col2
       fprintf('size different\n');
       return ;
   end
   [~, MseMinPos] = SearchByMSE(X1oriL, X2oriL) ;
   m = floor((MseMinPos(1) - 1) / 8) + 1;
   n = mod(MseMinPos(1) - 1, 8) + 1;
   % map  X1ori(5:row-12, 5:col-12) to X2ori(m:row-17+m, n:col-17+n)
   % find the links between 8x8 blocks
   X1L = X1oriL(5:rowL-4, 5:colL-4);
   X2L = X2oriL(m:rowL-9+m, n:colL-9+n) ;
   
%    [~, Cedges, blkcost_layer1] = DoEdmondMatch(X1L - mean(X1L(:)), X2L - mean(X2L(:)), 8) ;
   
%    ===========================================
   % in DCT
    T = DCT_X(8);
    
    C1 = ComputeDCTimage(X1L, T, 1) ;
    C2 = ComputeDCTimage(X2L, T, 1) ;
    
    mv = dlmread(mvfilename);
    [height, width] = size(C1(1 : 8 : end, 1 : 8 : end));
    idx = [1, 1; 1, 2; 2, 1; 3, 1; 2, 2; 1, 3];
    for i = 1
        py = reshape(mv(:, (i - 1) * 2 + 1), width, height)';
        px = reshape(mv(:, i * 2), width, height)';
        DCT1 = C1(idx(i, 1) : 8 : end, idx(i, 2) : 8 : end);
        DCT2 = C2(idx(i, 1) : 8 : end, idx(i, 2) : 8 : end);
        
        meanPx = mean(DCT1(:));
        D = abs(DCT1 - mean(DCT1(:)));
        hh = max(DCT1(:)) - meanPx;
        lh = min(DCT1(:)) - meanPx; 
%         mask = (DCT1 > (0.3 * hh + meanPx)) + (DCT1 < (0.3 * lh + meanPx));
%         mask = D > median(D(:));
        mask = DCT1 > median(DCT1(:));
%         ey = py - spatialPredict(py);
%         ex = px - spatialPredict(px);
%         EntropyArray((i - 1) * 4 + 1 : i * 4) = [calEntropy(py(:)), calEntropy(ey(:)), calEntropy(px(:)), calEntropy(ex(:))];
        close all;
        figure, 
        subplot(3, 2, 3)
        imshow(DCT1, []);
        subplot(3, 2, 4);
        imshow(DCT2, []);
        subplot(3, 2, 5);
        imshow(abs(py .* mask), [0, 26]);
%         imshow(mask, []);
        subplot(3, 2, 6);
        imshow(abs(px .* mask), [0, 36]);
        subplot(3, 2, 1);
        imshow(uint8(X1L));
        subplot(3, 2, 2);
        imshow(uint8(X2L));
    end

end

