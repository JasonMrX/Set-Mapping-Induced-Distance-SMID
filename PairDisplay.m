function PairDisplay(filename1, filename2, vote)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sizef = 264*264;
   %sizef = 136*136;
   [ X1oriL] = ReadResizeGrayImage(filename1, sizef) ;
   [ X2oriL] = ReadResizeGrayImage(filename2, sizef) ;
   %[row col] = size(X2ori);
   %X1ori(1:row, 1:col) = X1ori(1:row, 1:col)/0.95; % / max(max(X1ori(1:row, 1:col))) * max(max(X2ori(1:row, 1:col))) ;
   %X1ori(1:row/3, 1:col) = X2ori(1:row/3, 1:col)  ;
   %X2ori(1:row/3, 1:col) = X2ori(1:row/3, 1:col) -25 ;
   %X2ori(row/2:row, 1:col) = X2ori(row/2:row, 1:col) -15 ;
    
   [rowL, colL] = size(X1oriL);
   [row2, col2] = size(X2oriL);
   
   if rowL ~= row2 || colL ~= col2
        close all;
        figure,
        subplot(1, 2, 1);
        imshow(uint8(X1oriL));
        title(filename1)
        subplot(1, 2, 2);
        imshow(uint8(X2oriL));
        title(filename2)
   end

end

function [Gsmallcrop]= ReadResizeGrayImage(filename, sizef)
    I = imread(filename);
    
    n = size(I, 3);
    
    if n == 3
      G = rgb2gray(I);
    else
      G = I ;  
    end
    
    [row, col] = size(G);
     
    t = sqrt( sizef/row/col );
    Gsmall = double(imresize(G, t)) ;
    % crop it 
    [row, col] = size(Gsmall);
    Gsmallcrop = Gsmall(1:floor(row/8)*8, 1:floor(col/8)*8);
end