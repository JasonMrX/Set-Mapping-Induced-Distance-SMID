function [Y, MseMinPos] = SearchByMSE(X1ori, X2ori)
   [row, col] = size(X1ori);
   
   BlkMse = zeros(64,1);
   X1 = X1ori(5:row-4, 5:col-4);
   for m=1:8
   for n=1:8    
       X2 = X2ori(m:row-9+m, n:col-9+n);           
       BlkMse((m-1)*8+n) =  sum(sum((X1-X2).^2))*64/(row-8)/(col-8) ;
   end
   end
   [Y, MseMinPos] = sort(BlkMse);
end