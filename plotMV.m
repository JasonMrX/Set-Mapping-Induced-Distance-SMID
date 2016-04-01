argc = 5;
[ImagePair, OtherMetric] = ReadDatFile('ComputeResult_moreDat_P2x2_MVE_LB8x8_EMD_HI.dat', argc) ;

for i = 1 : numel(ImagePair)
   filename1 = ImagePair(i).name1;
   filename2 = ImagePair(i).name2;
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
       mDistDCT = 0;
       mDistDCTEntropy = 0;
       motionVectors = [];
       success = 0;
       return ;
   end
   success = 1;
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
    
    mvfilename = ['mv/', num2str(ImagePair(i).imageid)];
    a = dlmread(mvfilename);
    a1 = a(:, 3 : 4);
    plane1(:, :, 1) = reshape(a1(:, 1), 37, 27)';
    plane1(:, :, 2) = reshape(a1(:, 2), 37, 27)';
    py = plane1(:, :, 1);
    ey = py - spatialPredict(py);
    px = plane1(:, :, 2); 
    ex = px - spatialPredict(px);
    close all;
    figure, 
    subplot(2, 2, 1)
    imshow(C1(1 : 8 : end, 2 : 8 : end), []);
    subplot(2, 2, 2);
    imshow(C2(1 : 8 : end, 2 : 8 : end), []);
    subplot(2, 2, 3);
    imshow(py, [-26, 26]);
    subplot(2, 2, 4);
    imshow(px, [-36, 36]);
    
    title(num2str(ImagePair(i).imageid));
    
end


    