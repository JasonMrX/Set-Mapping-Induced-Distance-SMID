function DisplayFilesInClustersWithIDs(Nameclusterfilename, namelistfile)
% input is name of a cluster file, which shows filenames in clusters.
% the secodn input is a standard name list file name, which is used to
% generate the ID to be shown in the figures


%DisplayFilesInClustersWithIDs('AllNameCluster_revised_new0719.dat', 'AllFilelist_new0719.dat')
%Will generate images with clustered image together for verification.
%
%Should be run, whenever new clustering is used.    July-19, 2013

%save figure to images.

% this is for collecting clustering information. A IDcluster file is
% supposed to be generated. 

  namelist = ReadNamelistfile(namelistfile)   ;

  Cluster = ReadClusterfile(Nameclusterfilename) ;
  
  SaveToImage = 1 ;
  
  
    Font = GenerateFont ;

      ROW = 768*2 ;
      COL = 1024*2 ;
      Img(1:ROW, 1:COL, 1:3) = 0;

  
   N = length(Cluster);
   
   for n=1:N
       
      K = Cluster(n).size ;
      
      Img(1:ROW, 1:COL, 1:3) = 0;

      Kroot = ceil(sqrt(K));
      row = floor(ROW / Kroot) ;
      col = floor(COL / Kroot) ;
      
      
      for k=1:K
          
          
          filename = Cluster(n).name(k).str ;
          ID = findIDinlist(namelist, filename);
          
          I = imread(filename);
          I = imresize(I, [row, col]);
      
      
          idximg = PrintAnumber(Font, ID)*200;
          [idxrow idxcol]= size(idximg);
          I(2:idxrow+1, 2:idxcol+1) = ((double(I(2:idxrow+1, 2:idxcol+1)) + double(idximg))/2) ;
      
          Img(  floor((k-1)/Kroot)*row + 1 : floor((k-1)/Kroot)*row + row, mod((k-1), Kroot)*col + 1 : mod((k-1),Kroot)*col + col  , :  ) = I ;
      end

      if SaveToImage==1
          s = sprintf('Cluster_%d.jpg', n);
          imwrite(uint8(Img), s, 'Quality', 90);
      else
          figure
          imshow(uint8(Img));
          s = sprintf('Cluster %d', n);
          title(s);
      end
      
   end
      
      %imwrite(uint8(Img), strcat('ShowAll_', namelistfile, '.jpg'), 'Quality', 100) ;
      
function   ID = findIDinlist(namelist, filename)

  N = length(namelist);
  ID = 0 ;
  for n=1:N
      if strcmp(namelist(n).str, filename)
         ID = n-1;
         break ;
      end
  end
  
function Cluster = ReadClusterfile(filename)
fin = fopen(filename);
  
   if ~fin
       fprintf('Creating a new cluster file\n');
       
       fclose(fin);
       return ;
       
   end
   
    K = 0 ;
    count = 1;
    while count>0
      [M count] = fscanf(fin, '%d', 1) ;
      
      if count>0
          
        K = K + 1;
        
        Cluster(K).size = M ;        
        for m=1:M
          [name len] = fscanf(fin, '%s', 1) ;
          if len>0
              Cluster(K).name(m).str = name ;
          else
              fprintf('reading cluster file error\n');
          end
        end
        
      end
    end
fclose(fin);      
      
function namelist = ReadNamelistfile(namelistfile)  
    % read the name list
  fin = fopen(namelistfile,'r');
  N=1 ;
  len = 1;
  
  t = fscanf(fin, '%d', 1);
  
  while len>0
    [name len] = fscanf(fin, '%s ', 1) ;
    if len>0        
        namelist(N).str = name ;
        N = N + 1 ;
    end
  end
  fclose(fin);

  
function Font = GenerateFont

  n = 4 ;
  row = 5*n ;
  col = 4*n ;
  
  Font = zeros(row, col, 10);
  
  %for i=1:10
  %    Font(1:row, 1:col, i) = zeros(row, col);
  %end
  
  % 0
  digit = 0 ;
  
  
  Keys = [ 1 2 3 5 9  13 17 18 19 15 11 7];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 3 7 11 15 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [1 2 3 7 11 9 10 13 17 18 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 7 9 10 11 15 17 18 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 5 9 10 11 3 7 15 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 5 9 10 11 15 19 17 18];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 5 9 10 11 13 15 17 18 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 7 11 15 19 ];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 5 7 9 10 11 13 15 17 18 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  digit = digit + 1;
  
  Keys = [ 1 2 3 5 7 9 10 11 15 17 18 19];
  Font = DoOneDigit( digit, Keys, Font, n);
  
  
  
  
function FontOut = DoOneDigit( digit, keys, Font, n)
    FontOut = Font ;
    I = FontOut(:,:,digit+1) ;
    
    len = length(keys);
    
    for i=1:len
        BlockOut = FillOneSmallBlock(I, keys(i), n);
        I = BlockOut ;
    end
    
    FontOut(:,:,digit+1) = BlockOut;
    
    
  
  % size of Block is row - col
    function BlockOut = FillOneSmallBlock(Block, idx, n)
        
        %[row col] = size(Block);        
        BlockOut = Block ;
        BlockOut( floor(idx/4)*n+1:floor(idx/4)*n+n ,  mod(idx,4)*n+1: mod(idx,4)*n+n) = 1;
        
function Img = PrintAnumber(Font,  number ) 

 %Font = GenerateFont ;
 if number<0
     NumberOfDigit=0;
 end
 if number==0
     NumberOfDigit=1;
 else
     NumberOfDigit = 1 + floor( log10(number) ) ;
 end
 Img = zeros(24, NumberOfDigit*16+2);
 
 
 for i = 1: NumberOfDigit
    digitstr(i) =  mod(number,10);
    number = floor(number/10);
 end
 
 for i = 1: NumberOfDigit
     digit = digitstr(NumberOfDigit-i+1) ;
     Img(3:22,  (i-1)*16 +1:(i-1)*16 +16) = Font(:, :, digit+1) ;      
 end           
