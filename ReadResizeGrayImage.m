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
