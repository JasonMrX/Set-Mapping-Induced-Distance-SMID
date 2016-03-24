function  [mDistDCT, mDistDCTEntropy, motionVectors, success] = ComputeDCTComplexity_Entropy(filename1, filename2)
% To compute ISmai with two-layer blk matching.
% First layer uses non-overlapping perfect match
% Second layer uses co-located area, but allow smaller blks.

% modified based on 7, to see what we get with larger image size without
% sliding window.

% within 8x8 block, try to match pixel by pixel

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
    
%     ================= mv entropy ====================
    mv = [];
    motionVectors = [];
    idx = [1, 1; 1, 2; 2, 1; 3, 1; 2, 2; 1, 3];
    
    for i = 1 : size(idx, 1)
        DCTImage1 = C1(idx(i, 1) : 8 : end, idx(i, 2) : 8 : end);
        DCTImage2 = C2(idx(i, 1) : 8 : end, idx(i, 2) : 8 : end);
        [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
        newMVs = getMVFromEdges(Cedges, size(DCTImage1, 2));
        mv = [mv; newMVs];
        motionVectors = [motionVectors newMVs];
    end
    
%     DCTImage1 = C1(1 : 8 : end, 1 : 8 : end);
%     DCTImage2 = C2(1 : 8 : end, 1 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors = repmat(mv, 1, 6);
%     
%     DCTImage1 = C1(1 : 8 : end, 2 : 8 : end);
%     DCTImage2 = C2(1 : 8 : end, 2 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors(:, 3 : 4);
%     
%     DCTImage1 = C1(2 : 8 : end, 1 : 8 : end);
%     DCTImage2 = C2(2 : 8 : end, 1 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors(:, 5 : 6);
%     
%     DCTImage1 = C1(3 : 8 : end, 1 : 8 : end);
%     DCTImage2 = C2(3 : 8 : end, 1 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors(:, 7 : 8);
%     
%     DCTImage1 = C1(2 : 8 : end, 2 : 8 : end);
%     DCTImage2 = C2(2 : 8 : end, 2 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors(:, 9 : 10);
%     
%     DCTImage1 = C1(1 : 8 : end, 3 : 8 : end);
%     DCTImage2 = C2(1 : 8 : end, 3 : 8 : end);
%     [~, Cedges, ~] = DoEdmondMatch(DCTImage1 - mean(DCTImage1(:)), DCTImage2 - mean(DCTImage2(:)), 1) ;
%     mv = [mv; getMVFromEdges(Cedges, size(DCTImage1, 2))];
%     motionVectors(:, 11 : 12);
    
    mDistDCTEntropy = calEntropy(mv(:));
%     ================= mv entropy ====================

    mDistDCT = ComputeEMD4TwoDCTImage(C1, C2)/64 ;
    
function [mv] = getMVFromEdges(Cedges, col)
    Id1 = Cedges(:, 1);
    Id2 = Cedges(:, 2);
    Pos1 = [floor(Id1 / col) + 1, mod(Id1, col) + 1];
    Pos2 = [floor(Id2 / col) + 1, mod(Id2, col) + 1];
    mv = Pos1 - Pos2;
    
function [entropy] = calEntropy(data)
    h = hist(data, min(data):max(data));
    h = h(h ~= 0);
    h = h / sum(h);
    entropy = sum(-h .* log2(h));
      
% search by MSE for a few locations 
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
    

function [blkcost, distcost, motionVectors] = DoEmdondMatchLayer2(edges, X1, X2)
    Ne = size(edges, 1);
    [row, col] = size(X1);

    motionVectors = zeros(numel(X1) / 4, 2);
    blkrowL = row / 8;
    blkcolL = col / 8;   

    totalcost = 0;
    distcost = 0;
    
    for n=1:Ne
        Id1 = edges(n,1);
        Id2 = edges(n,2);
        
        blki1 = floor(Id1/blkcolL) + 1 ;
        blkj1 = mod(Id1, blkcolL) +1 ;
        blki2 = floor(Id2/blkcolL) + 1 ;
        blkj2 = mod(Id2, blkcolL) +1 ;
        
        Center_i1 = (blki1-1)*8 + 1 ;
        Center_j1 = (blkj1-1)*8 + 1 ;
        Center_i2 = (blki2-1)*8 + 1 ;
        Center_j2 = (blkj2-1)*8 + 1 ;
        
        block1 = X1( Center_i1:Center_i1+7, Center_j1:Center_j1+7 );
        block2 = X2( Center_i2:Center_i2+7, Center_j2:Center_j2+7 );
        [mv, cost] = DoHungarianMatch(block1, block2, 2) ;
        Pos1 = [floor(Id1 / blkcolL) + 1, mod(Id1, blkcolL) + 1];
        Pos2 = [floor(Id2 / blkcolL) + 1, mod(Id2, blkcolL) + 1];
        mv = repmat(Pos2 - Pos1, size(mv, 1), 1) * 8 + mv;
        totalcost = totalcost + cost * 16 ;
        distcost = distcost + 4 * sum(sqrt(mv(:, 1).^2 + mv(:, 2).^2));
        
        motionVectors((n - 1) * 16 + 1 : n * 16, :) = mv;
    end
    
    blkcost = totalcost / Ne ;  % per block
    distcost = distcost / Ne;
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ClinkMat, Clink, cost] = DoEdmondMatch(X1, X2, blksize)
      [row, col] = size(X1);
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

      scale = 100 ;
      % from weight matrxi to graph
      weight = round(weight*scale) ;  % convert to integers
      Graph = zeros(Len*Len, 3) ;
      pos = 1; 
      for m=1:Len
      for n=1:Len
          Graph(pos, 1:3) = [ m-1, n-1+Len,  weight(m,n) ] ;
          pos = pos + 1;
      end
      end

    WriteGraphBindata('tmpdat.bin', Len, Len*Len, Graph) ;

    system('blossom.exe -e tmpdat.bin -w tmplink.bin');


    [~, Clink, ClinkMat, matchcost] = ReadLinks('tmplink.bin') ;

  
      %[link cost] = Hungarian(weight);
      
      cost = matchcost/Len/scale;  % per block 
      
      
      if 0
                  if Len<=256*256
                    [link Huncost] = Hungarian(weight);
                    if Huncost ~= matchcost
                        fprintf('Hungarian does not match with Edmond. Something wrong.\n');
                    else
                        fprintf('Hungarian matches with Edmond.\n');
                    end
                  end
      end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mv, cost] = DoHungarianMatch(X1, X2, blksize)
      [row1, col1] = size(X1);
      [row2, col2] = size(X2);
      row = min(row1, row2);
      col = min(col1, col2);
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
      [link, cost] = Hungarian(weight);
      cost = cost/Len;  % per block 
      mv = link2mv(link, blksize);

%%%%%%%%%%%%%%%%%%%%% transform link matrix to motion vection mv %%%%%%%%%%
function [mv] = link2mv(link, blksize)
    numVertex = size(link, 1);
    blkWidth = sqrt(numVertex);
    idx1 = mod(find(link(:)) - 1, numVertex) + 1;
    idx2 = (1 : numVertex)';
    pos1 = [floor((idx1 - 1) / blkWidth) + 1, mod(idx1 - 1, blkWidth) + 1];
    pos2 = [floor((idx2 - 1) / blkWidth) + 1, mod(idx2 - 1, blkWidth) + 1];
    mv = (pos2 - pos1) * blksize;
    
 
%%%%%%%%%%%%%
% Clink stores (id1 id2)-pairs of edges, idx\in[0, N-1]
function   [Nnode, Clink, ClinkMat, matchcost] = ReadLinks(filename)
  fin = fopen(filename, 'r');
  
  Nnode = fscanf(fin, '%d', 1);
  Nlink = fscanf(fin, '%d', 1);
  
  Clink = zeros(Nlink, 2);
  tmp = fscanf(fin, '%d', 2*Nlink);
  for n=1:Nlink
      Clink(n,1) = tmp( n*2 -1 ) ;
      Clink(n,2) = tmp( n*2 ) ;
  end
  fclose(fin);
  
  
  % convert to matrix format for matlab use.
  ClinkMat=zeros(Nnode/2, Nnode/2);
  SingleNnode = Nnode/2 ;
  for n=1:Nlink % in C code, diff indexes are used for col and row nodes.
     Clink(n,2) = Clink(n,2) - SingleNnode ; % In C, we need to differ, but not here
     ClinkMat( Clink(n,1)+1, Clink(n,2)+1 ) = 1 ; 
  end

  fin=fopen('matchcost.dat', 'r');
  matchcost = fscanf(fin, '%f', 1);
  fclose(fin);

  
% GraphEdges store Nedge-x-3 matrix for the edges.
function WriteGraphBindata(filename, Nnode, Nedge, GraphEdges)
   fout = fopen(filename, 'wb');
   fwrite(fout, 2*Nnode, 'integer*4');
   fwrite(fout, Nedge, 'integer*4');
   Graphout = GraphEdges' ;  % write column by column
   fwrite(fout, Graphout, 'integer*4');
   fclose(fout);  

%    ==================== DCT funcs ========================
function T = DCT_X(Dim)
    N=Dim;
    M=Dim;
    T = zeros(Dim);
    for u=0:Dim-1
        for i=0:Dim-1
            T(u+1,i+1) = cos(pi*u*(2*i+1)/2/N) * Gam(u) * sqrt(2/N) ;
        end
    end
    
function x = Gam(i)
    if i==0
        x = 1/sqrt(2);
    else
        x = 1 ;
    end
    
function C = ComputeDCTimage(X, T, type)
% type = 1 : forward
%      = 2 : inverse
    [row, col] = size(X) ;
    N = size(T, 1);
    if type == 2
        T = T';
    end
    for j = 1 : N : row
        for i = 1 : N : col
            block = X(j : j + N - 1, i : i + N - 1) ;
            C(j : j + N - 1, i : i + N - 1) = (T * block) * (T') ;
        end
    end
    
function TotalEmdDCT = ComputeEMD4TwoDCTImage(C1, C2)    
    TotalEmdDCT = 0 ;
    [row1, col1] = size(C1);
    [row2, col2] = size(C2);
    for i = 1 : 8
        for j = 1 : 8
            cblock1 = C1(i : 8 : row1, j : 8 : col1);
            cblock2 = C2(i : 8 : row2, j : 8 : col2);

            Seq1 = cblock1(:);
            Seq2 = cblock2(:);

            [emd, emdNew, ~] = computeEMD4TwoSeq(Seq1, Seq2, 0) ; %#ok<ASGLU>

            %if i~=1 || j~=1
            TotalEmdDCT = TotalEmdDCT + emdNew;    
            %end
            %fprintf('%6.1f ', emdNew);
        end
    end

function [emd, emdNew, vard] = computeEMD4TwoSeq(Seq1, Seq2, ShowFigure)
    sample = min(min(Seq1), min(Seq2)) : 1 : max( max(Seq1), max(Seq2) ) ;
    p1 = hist(Seq1, sample);
    p1 = p1/sum(p1);
    p2 = hist(Seq2, sample);
    p2 = p2/sum(p2);

    %ShowFigure=0; 
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

    vard = sum(abs(p1 - p2));   
    %fprintf('Distance used to be %s %s   emd=%f     varDis=%f    %f\n', filename1, filename2, emd, vard,  sum(abs(cdf1- cdf2)));

    Len = length(p1);
    E1 = sum((1 : Len) .* p1) ;
    E2 = sum((1 : Len) .* p2) ;

    %fprintf('E1=%f E2=%f\n', E1, E2);
    if E1 < E2
        shift = round(E2 - E1) ;
        p1new(1 : Len + shift)=0;
        p2new(1 : Len + shift)=0;
        p2new(1 : Len) = p2;
        p1new(shift + 1 : Len + shift) = p1 ;
        p1 = p1new;
        p2 = p2new;
    else
        shift = round(E1 - E2) ;
        p1new(1 : Len + shift)=0;
        p2new(1 : Len + shift)=0;
        p1new(1 : Len)=p1;
        p2new(shift + 1 : Len + shift) = p2 ;
        p1 = p1new;
        p2 = p2new;
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
    c(1 : N) = cdf1;
    c(N + 1 : 2 * N) = cdf2 ;
    c = sort(c);
    
    Thr(1,1 : N) = cdf1;
    Thr(2,1 : N) = cdf2;
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
        if c(n) >= ThMax(posmax) && (posmax<N)
            L = L+1;
            posmax = posmax + 1;
        end
        if c(n) >= ThMin(posmin) && (posmin<N)
            R = R + 1;
            posmin = posmin + 1;
        end
    end
    
function cdf = computeCdfFromHist(p)
    p = p / sum(p); % in case inputs are histogram
    Nref = length(p);
    cdf = zeros(size(p));
    cdf(1) = p(1);
    for n = 2 : Nref
        cdf(n) = cdf(n - 1) + p(n);
    end