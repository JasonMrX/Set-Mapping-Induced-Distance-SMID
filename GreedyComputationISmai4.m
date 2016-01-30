function  [success hist_IS mDistPixel mDistDCT pmse] = GreedyComputationISmai4(filename1, filename2)
% Downsample to certain size, find 64 sliding images and match

   
   sizef = 136*136;
   X1ori = ReadResizeGrayImage(filename1, sizef) ;
   X2ori = ReadResizeGrayImage(filename2, sizef) ;
   
   % watch the same size
   [row col] = size(X1ori);
   [row2 col2] = size(X2ori);
   if row~=row2 || col~=col2
       fprintf('size different\n');
       success = 0 ;
       hist_IS = 0 ;
       mDistPixel=0;
       mDistDCT=0;
       pmse=0;
       return ;
   end
   
   success = 1 ;
 
   % histogram intersection
   hist_IS = ComputeHistogramIntersection(X1ori(:), X2ori(:));
   
   T = DCT_X(8);
   
   [BlkMse MseMinPos] = SearchByMSE(X1ori, X2ori) ;
   
   K = 6 ;
   DistPixel = zeros(1,K);
   DistDCT   = zeros(1,K);
   
   X1 = X1ori(5:row-4, 5:col-4);
   C1 = ComputeDCTimage(X1, T, 1) ;
   for k=1:6
       m = floor((MseMinPos(k)-1)/8)+1 ;
       n = mod( MseMinPos(k)-1, 8) + 1;
       X2 = X2ori(m:row-9+m, n:col-9+n);     
       DistPixel(k) = DoHungarianMatch(X1-mean(mean(X1)), X2-mean(mean(X2)), 8) ;
       %DistPixel(m,n) = DoHungarianMatch(X1, X2, 8) ;
       
       %fprintf('meanX1=%8.2f   meanX2=%8.2f    \n', mean(mean(X1)), mean(mean(X2)));
       
       C2 = ComputeDCTimage(X2, T, 1) ;
       
       DistDCT(k) = ComputeEMD4TwoDCTImage(C1, C2) ;
       %DistDCT(m,n)=0;
       
       fprintf('%d %d    %8.2f   %8.2f       Mse=%8.2f\n', m,n, DistPixel(k)/64,  DistDCT(k)/64, BlkMse(k)/64  );
   end

   mDistPixel = (min(DistPixel))/64;
   mDistDCT = (min(DistDCT))/64;
   pmse = (min(BlkMse))/64;
   
   fprintf('PixelBlkMatch Dist=%8.2f     DCTBlkMatch=%8.2f  (max=%8.2f)    mse=%8.2f\n', (min(DistPixel))/64, (min(DistDCT))/64, (max(DistDCT))/64, (min(BlkMse))/64);
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   hist_IS = ComputeHistogramIntersection(seq1, seq2)
    miny = min(min(seq1), min(seq2));
    maxy = max(max(seq1), max(seq2));
    Y1 = hist(seq1, miny:maxy);
    Y2 = hist(seq2, miny:maxy);
   
    p1 = Y1/sum(Y1);
    p2 = Y2/sum(Y2);
    
    %hist_IS = sum(min(p1, p2));
    
    vard = sum(abs(p1-p2));
    
    hist_IS = 1-vard/2 ;
    
    %fprintf('%f   %f     %f\n', hist_IS, vard,   1-vard/2);
    
    
% search by MSE for a few locations 
function [Y MseMinPos] = SearchByMSE(X1ori, X2ori)
   [row col] = size(X1ori);
   
   BlkMse = zeros(64,1);
   X1 = X1ori(5:row-4, 5:col-4);
   for m=1:8
   for n=1:8    
       X2 = X2ori(m:row-9+m, n:col-9+n);           
       BlkMse((m-1)*8+n) =  sum(sum((X1-X2).^2))*64/(row-8)/(col-8) ;
   end
   end
   [Y, MseMinPos] = sort(BlkMse);     
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function TotalEmdDCT = ComputeEMD4TwoDCTImage(C1, C2)    
     %%%%%%%%%%%%%%%%%%%%%%%%Compute in DCT domain%%%%%%%%%%%%%%%%%%%%%%%
      TotalEmdDCT = 0 ;
      [row1 col1] = size(C1);
      [row2 col2] = size(C2);
      for i=1:8
      for j=1:8

          cblock1 = C1(i:8:row1, j:8:col1);
          cblock2 = C2(i:8:row2, j:8:col2);

          Seq1 = cblock1(:);
          Seq2 = cblock2(:);

          [emd emdNew vard] = computeEMD4TwoSeq(Seq1, Seq2) ;
          TotalEmdDCT = TotalEmdDCT + emdNew ;    
          %fprintf('%6.1f ', emdNew);
      end
      end
      
      %fprintf('  EmdDCT=%f    mse=%f\n',  TotalEmdDCT, sum(sum((L1_X1 - L1_X2).^2))/(row*col/64));
      
function Gsmallcrop = ReadResizeGrayImage(filename, sizef)
    I = imread(filename);
    
    [r c n] = size(I);
    
    if n==3
      G = rgb2gray(I);
    else
      G = I ;  
    end
    
    [row col] = size(G);
     
    t = sqrt( sizef/row/col );
    Gsmall = double(imresize(G, t)) ;
    % crop it 
    [row col] = size(Gsmall);
    Gsmallcrop = Gsmall(1:floor(row/8)*8, 1:floor(col/8)*8);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost] = DoHungarianMatch(X1, X2, blksize)
      [row col] = size(X1);
      M = round(row/blksize);
      N = round(col/blksize);
      
      Len = M*N;
         
      weight = zeros(Len, Len);
      for r1=1:M
      for c1=1:N
          block1 = X1( (r1-1)*blksize+1:r1*blksize, (c1-1)*blksize+1:c1*blksize );
          for r2=1:M
          for c2=1:N
             block2 = X2( (r2-1)*blksize+1:r2*blksize, (c2-1)*blksize+1:c2*blksize );
             weight((r1-1)*N+c1, (r2-1)*N+c2) = sum(sum(((block1-block2).^2)));
          end
          end
      end
      end
      [link cost] = Hungarian(weight);
      
      cost = cost/Len;  % per block 
      
function blockP = FindAblockByPred(block, X, h, w, N, Range)
    % search in a neighborhood for a best match of block in X
    %Range=4;
    [row col] = size(X);
    dmin = 2^30;
    them=h;
    then=w;
    for m=max(h-Range,1):min(h+Range, row-N+1)
    for n=max(w-Range,1):min(w+Range, col-N+1)
        Pblock = X(m:m+N-1, n:n+N-1);
        d = sum(sum((block-Pblock).^2)) ;
        if d<dmin
            dmin = d;
            them = m;
            then = n;
        end
    end
    end
    blockP = X(them:them+N-1, then:then+N-1);
        
    

% X1 X2 are pre-processed to have the same size.
function [stat1 stat2, totalD] = SmallBlockMatchCost(X1, X2, BlkWidth)
  
  [row col] = size(X1);
  
  newX1 = X1 ;
  newX2 = X2 ;      
  
  blkRow = row/BlkWidth ;
  blkCol = col/BlkWidth ;
  v1 = zeros(blkRow*blkCol,1 );
  v2 = zeros(blkRow*blkCol,1 );
  for k = 1:blkRow*blkCol
      %m = floor((k-1)/blkCol) + 1 ;
      %n = mod(k-1, blkCol) + 1 ;
      %im = (m-1)*BlkWidth +1 ;
      %in = (n-1)*BlkWidth +1 ;
      %blkX1 = newX1(im:im+BlkWidth-1, in:in+BlkWidth-1);
      %blkX2 = newX2(im:im+BlkWidth-1, in:in+BlkWidth-1);   
      blkX1 = newX1( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );
      blkX2 = newX2( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );
      v1(k) = sum(sum(blkX1.*blkX1));
      v2(k) = sum(sum(blkX2.*blkX2));
  end

  TotalV(1:blkRow*blkCol)=v1 ;
  TotalV(blkRow*blkCol+1: blkRow*blkCol*2)=v2 ;
  
  [En, Id]= sort(TotalV, 'descend');

  stat1(1:blkRow*blkCol) = -1 ;  % -1 means available, otherwise would be >=1, as index
  stat2(1:blkRow*blkCol) = -1 ;
  totalD = 0 ;
  for n=1:length(Id)
     idx = Id(n);     
     if idx>blkRow*blkCol
        k = idx - blkRow*blkCol ;
        if stat2(k)>0  % matched already
            continue;
        end
        blk = newX2( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );

        [d theK] = FindAmatch(blk, newX1, blkRow, blkCol, BlkWidth, stat1);
        stat1(theK) = k ;
        stat2(k) = theK ;
        totalD = totalD + d ;
     else
        k = idx ; 
        if stat1( k)>0
            continue;
        end 
        blk = newX1( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );

        [d theK] = FindAmatch(blk, newX2, blkRow, blkCol, BlkWidth, stat2);
        stat1(k) = theK ;
        stat2(theK) = k ;
        totalD = totalD + d ;
     end     
      
  end
  
  % verfiy the calculation
  totalDby1 = 0 ;
  for k = 1:blkRow*blkCol
      theK = stat1(k) ;
      
      blkX1 = newX1( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );
      blkX2 = newX2( floor((theK-1)/blkCol)*BlkWidth+1:(floor((theK-1)/blkCol)+1)*BlkWidth, mod(theK-1,blkCol)*BlkWidth+1:(mod(theK-1,blkCol)+1)*BlkWidth );
      d = sum(sum((blkX1-blkX2).^2));
      totalDby1 = totalDby1 +d ; ;
  end

  totalDby2 = 0 ;
  for k = 1:blkRow*blkCol
      theK = stat2(k) ;
      
      blkX1 = newX2( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );
      blkX2 = newX1( floor((theK-1)/blkCol)*BlkWidth+1:(floor((theK-1)/blkCol)+1)*BlkWidth, mod(theK-1,blkCol)*BlkWidth+1:(mod(theK-1,blkCol)+1)*BlkWidth );
      d = sum(sum((blkX1-blkX2).^2));
      totalDby2 = totalDby2 +d ; ;
  end
  
  if round(totalDby1)~=round(totalDby2) || round(totalDby1)~=round(totalD)
      fprintf('Dby1=%f   Dby2=%f   totalD=%f\n', totalDby1, totalDby2, totalD);
  end
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d theK] = FindAmatch(blk, newX, blkRow, blkCol, BlkWidth, stat)
  minD = 2^30 ;
  mink = -1 ;
  for k = 1:blkRow*blkCol
      %m = floor((k-1)/blkCol) + 1 ;
      %n = mod(k-1, blkCol) + 1 ;
      %im = (m-1)*BlkWidth +1 ;
      %in = (n-1)*BlkWidth +1 ;
      %blkX1 = newX1(im:im+BlkWidth-1, in:in+BlkWidth-1);
      %blkX2 = newX2(im:im+BlkWidth-1, in:in+BlkWidth-1);   
      if stat(k)<0 
        blkX = newX( floor((k-1)/blkCol)*BlkWidth+1:(floor((k-1)/blkCol)+1)*BlkWidth, mod(k-1,blkCol)*BlkWidth+1:(mod(k-1,blkCol)+1)*BlkWidth );
        d = sum(sum( (blkX - blk).^2)) ;
        if d< minD
            minD = d ;
            mink = k ;
        end
      end      
  end
  theK = mink ;
  d = minD ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [emd emdNew vard] = computeEMD4TwoSeq(Seq1, Seq2)
   %sample = -2048:1:2048 ;
   sample = min(min(Seq1), min(Seq2) ) : 1 : max( max(Seq1), max(Seq2) ) ;
   %fprintf('(%d %d %d %d)\n', max(Seq1), max(Seq2), min(Seq1), min(Seq2));
   %sample = -4096:1:4096 ;
   p1 = hist(Seq1, sample);
   p1 = p1/sum(p1);
   p2 = hist(Seq2, sample);
   p2 = p2/sum(p2);
   
   ShowFigure=0; 
       
   if ShowFigure==1 
   figure
   plot(sample, p1, 'r.-');
   hold on
   plot(sample, p2, 'k.-');
   %hold on 
   legend('filename1', 'filename2');
   end
   
   rou = 2 ;
           cdf1 = computeCdfFromHist(p1);
           cdf2 = computeCdfFromHist(p2);      
           emd = ComputeEMD(cdf1, cdf2, rou);   
           
           vard = sum(abs(p1-p2));   
           %fprintf('Distance used to be %s %s   emd=%f     varDis=%f    %f\n', filename1, filename2, emd, vard,  sum(abs(cdf1- cdf2)));

   
   Len=length(p1);
   E1 = sum((1:Len).*p1) ;
   E2 = sum((1:Len).*p2) ;
   
   %fprintf('E1=%f E2=%f\n', E1, E2);
   
   
   if E1<E2
       shift = round(E2-E1) ;
       p1new(1:Len+shift)=0;
       p2new(1:Len+shift)=0;
       p2new(1:Len)=p2;
       p1new(shift+1 :Len+shift) = p1 ;
       
       p1=p1new;
       p2=p2new;
   else
       shift = round(E1-E2) ;
       p1new(1:Len+shift)=0;
       p2new(1:Len+shift)=0;
       p1new(1:Len)=p1;
       p2new(shift+1 :Len+shift) = p2 ;
       
       p1=p1new;
       p2=p2new;
   end
   
   
   cdf1 = computeCdfFromHist(p1);
   cdf2 = computeCdfFromHist(p2);
      
   
  if ShowFigure==1 
  len=length(cdf1);
   sample = 1:len ;   
   figure
   stem(sample, cdf1, 'r.-');
   hold on
   stem(sample, cdf2, 'k.-');
   hold on 
   legend('filename1', 'filename2');
  end
  
   emdNew = ComputeEMD(cdf1, cdf2, rou);   

       
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
  
function cdf = computeCdfFromHist(p)
  p = p/sum(p); % in case inputs are histogram
  Nref = length(p);
  cdf(1) = p(1);
  for n=2:Nref
      cdf(n) = cdf(n-1) + p(n);
  end
    
function [newC] = GetDCTimage(filename, Dim)                  
  %I = imread(strcat(filename,'.jpg'));
  %figure
  %imshow(I);
    command = sprintf('jpeg_dec.exe -outfile x.img %s', strcat(filename,'.jpg'));
      
    system(command);
    % get C   
    [Q  Ntable] = readTQ ;
    
    [uY uU uV] = readU(Ntable) ;

      Chat_y = DequanImage(uY,Q(:,:,1)) ;
      Chat_u = DequanImage(uU,Q(:,:,2)) ;
      Chat_v = DequanImage(uV,Q(:,:,2)) ;

    
    Chat_y = AdjustDC(Chat_y) ;
    Chat_u = AdjustDC(Chat_u) ;
    Chat_v = AdjustDC(Chat_v) ;
    
    if Dim==8
        newC = Chat_y;
    else
        T = DCT_X(8);  
        X = ComputeDCTimage(Chat_y, T, 2) ;
        newXf = round(imresize(X, 8/Dim));
        [row col] = size(newXf);
        newX = newXf(1:floor(row/Dim)*Dim, 1:floor(col/Dim)*Dim) ;
        newC = round(ComputeDCTimage(newX, T, 1)) ;
    end
    %figure
    %imshow(uint8(DCx.*(DCx<50)));
    
    %figure
    %[y , I] = hist(DCx, 0:255);
    %plot(y/sum(y), 'r.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function It= AdjustSizeTo8plex(I)
[row col noc] = size(I);
    
newrow = floor(row/16)*16;
newcol = floor(col/16)*16;

It = I(1:newrow, 1:newcol, :);



function OLyLP = lowpass(U,a,b)
    [row, col] = size(U);
    LPMcut = [1 1 1 1 0 0 0 0; 1 1 1 0 0 0 0 0; 1 1 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;];
    
    fprintf('a=%f\n',a);
    LPM = lowpassMatrix(8, a,b).*LPMcut
    for i=1:8:row
    for j=1:8:col
        OLyLP(i:i+7, j:j+7) = U(i:i+7, j:j+7) .* LPM;
    end
    end
    
function OLyLP = lowpass22(U)
    [row, col] = size(U);    
    LPM = [1 1 1 1 0 0 0 0; 1 1 1 0 0 0 0 0; 1 1 0 0 0 0 0 0; 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 ;];
    for i=1:8:row
    for j=1:8:col
        OLyLP(i:i+7, j:j+7) = U(i:i+7, j:j+7) .* LPM;
    end
    end

    
function LPM = lowpassMatrix(N, a,b)    
  for i=0:N-1
  for j=0:N-1
    LPM(i+1,j+1) = (a+2*cos(i*pi/8))*(b+2*cos(j*pi/8))/(a+2)/(b+2);
  end
  end
    
function    DisOLsign_block = GetWholeOLblock(DisOLsign)

    thr = 0;
    [row col] = size(DisOLsign) ;
    
    DisOLsign_block = zeros(row, col);
    
    for i=1:8:row
    for j=1:8:col
        block = DisOLsign(i:i+7,j:j+7);
        
        if sum(block(:))>thr
            DisOLsign_block(i:i+7, j:j+7) = ones(8,8);
        end
    end
    end
    
    
function [uY uU uV] = readU(Ntable)
  finputU    = fopen('deJPEG.U','r');

  My = fread(finputU, [1, 1], 'uint32');
  Ny = fread(finputU, [1, 1], 'uint32');

  if Ntable>1
      Mu = fread(finputU, [1, 1], 'uint32');
      Nu = fread(finputU, [1, 1], 'uint32');

      Mv = fread(finputU, [1, 1], 'uint32');
      Nv = fread(finputU, [1, 1], 'uint32');
  end
  
  Uy = fread (finputU, [My, Ny], 'int16') ;

  if Ntable>1
    Uu = fread (finputU, [Mu, Nu], 'int16') ;
    Uv = fread (finputU, [Mv, Nv], 'int16') ;
  end
  
  fclose(finputU);
  uY = Uy' ;
  if Ntable>1
    uU = Uu' ;
    uV = Uv' ;
  else
    uU = 0;
    uV = 0;
  end
  
function [Q Ntable]= readTQ
  finput    = fopen('deJPEG.TQ','r');

  Len = fread(finput, [1,1], 'uint8');
  z = fread(finput, [1, Len], 'uint8') ;
  Ntable = z(1+z(1)*2+1) ;
  for n=1:Ntable
      Qtmp(1, 1:64) = z(z(1)*2+3+(n-1)*64: z(1)*2+2+n*64) ;
      Q(1:8,1:8,n) = [Qtmp(1:8);Qtmp(9:16);Qtmp(17:24);Qtmp(25:32);Qtmp(33:40);Qtmp(41:48);Qtmp(49:56);Qtmp(57:64)] ;
  end

  fclose(finput);
          
function Chat = DequanImage(U,Q)
[row, col] = size(U) ;
[N, N] = size(Q);


  for j = 1:N:row
      for i = 1:N:col
          block = U(j:j+N-1, i:i+N-1) ;
          %T
          %(T * block) * (T')
          Chat(j:j+N-1, i:i+N-1) = (Q .* block) ;
      end
  end    
function C = ComputeDCTimage(X, T, type)
% type = 1 : forward
%      = 2 : inverse
[row, col] = size(X) ;
[N, N] = size(T);

if type==2
    T = T' ;
end

  for j = 1:N:row
      for i = 1:N:col
          block = X(j:j+N-1, i:i+N-1) ;
          %T
          %(T * block) * (T')
          C(j:j+N-1, i:i+N-1) = (T * block) * (T') ;
      end
  end
  
function T = DCT_X(Dim)

  N=Dim;
  M=Dim;
  
  for u=0:Dim-1
  for i=0:Dim-1
    T(u+1,i+1) = cos(pi*u*(2*i+1)/2/N) * Gam(u) * sqrt(2/N) ;
  end
  end
  
  
  
  for j=0:Dim-1
  for v=0:Dim-1
    V(j+1,v+1) = cos(pi*v*(2*j+1)/2/M) * Gam(v) * sqrt(2/M) ;
  end
  end
  
  %T-V'
  %T
    
  %T*(T')
 
 function x = Gam(i)
  if i==0
    x = 1/sqrt(2);
  else
    x = 1 ;
  end
function  C_adjusted = AdjustDC(C) 
    [row col] = size(C);
    C_adjusted = C ;
    C_adjusted(1:8:row, 1:8:col) = C(1:8:row, 1:8:col) + 128*8 ;  