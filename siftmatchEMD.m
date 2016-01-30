function [matches score]= siftmatchEMD(descrip1, descrip2, match_thr)

  % make them a real histogram
  descriptH1 = EnsureHistogram(descrip1);
  descriptH2 = EnsureHistogram(descrip2);
  
  cdf1 = computeCdfFromPdf(descriptH1);
  cdf2 = computeCdfFromPdf(descriptH2);
  
  [t  M] = size(cdf1);
  [t  N] = size(cdf2);
 
  if M==0 || N==0
      score(1) = 2 ; % a large value ;
      matches(1,1) = 1 ;
      matches(2,1) = 1 ;
      return ;
  end
  
  Dist=zeros(M,N);
  
  for i=1:M
  for j=1:N
      Dist(i,j) =  sum( abs(cdf1(:,i) - cdf2(:,j))) ; % EMD_L1
      %Dist(i,j) = sum((descrip1(:,i) - descrip2(:,j)).^2) ;
  end
  end
  
  thre = match_thr ;
  K = 1 ;
          % Match from I2 to I1
         % [V, pos] = min(Dist, [], 1);
         % Disttmp = Dist ;
         % for i=1:N
         %    Disttmp(pos(i),i) = 10000 ;
         % end
         % [Vsec, possec] = min(Disttmp, [], 1);
         % Vsec = Vsec / thre ;
         %
         % for i=1:N
         %     if V(i)<Vsec(i)
         %         matches(1, K) = pos(i);
         %         matches(2, K) = i;
         %         score(K) = V(i);
         %         K = K + 1;
         %     end
         % end
  
  %Match from I1 to I2
  [V, pos] = min(Dist, [], 2);
  Disttmp = Dist ;
  for i=1:M
     Disttmp(i, pos(i)) = 10000 ;
  end
  [Vsec, possec] = min(Disttmp, [], 2);
  Vsec = Vsec / thre ;
  
  for i=1:M
      if V(i)<Vsec(i)
          matches(1,K) = i;
          matches(2,K) = pos(i);
          score(K) = V(i);
          K = K + 1;
      end
  end
  
  if K==1
      score(1) = 2 ; % a large value ;
      matches(1,1) = 1 ;
      matches(2,1) = 1 ;
  end
  
function descriptH = EnsureHistogram(descript)
  [row col] = size(descript);
  S = sum(descript, 1);
  descriptH = zeros(row, col);
  for c = 1:col
      descriptH(:,c) = descript(:,c)/S(c);
  end
  
  
function cdf = computeCdfFromPdf(descriptH)
  [row col] = size(descriptH) ;
  
  cdf(1,:) = descriptH(1, :);
  for n=2:row
      cdf(n,:) = cdf(n-1,:) + descriptH(n,:);
  end
  
function emd = ComputeEMD_L1(cdf1, cdf2)
  emd = sum(abs(cdf1-cdf2));

function emd = ComputeEMD(cdf1, cdf2, rou)
    N = length(cdf1);
    c(1:N) = cdf1;
    c(N+1:2*N) = cdf2 ;
    c = sort(c);
    
    Thr(1,1:N) = cdf1;
    Thr(2,1:N) = cdf2;
    ThMax = max(Thr, [], 1);
    ThMin = min(Thr, [], 1);
   
    L = 0 ;
    R = 1 ;
    emd = 0 ;
    %rou = 0.5;
    posmin = 2 ;
    posmax = 1 ;
    
    for n=2:2*N-2
        emd = emd + (c(n)-c(n-1)) * ((R-L)^rou) ;
        if c(n)>= ThMax(posmax) && (posmax<N)
            L = L+1;
            posmax = posmax + 1;
        end
        if c(n)>= ThMin(posmin)  && (posmin<N)
            R = R+1;
            posmin = posmin + 1;
        end
    end  
  