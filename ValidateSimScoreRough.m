function ValidateSimScoreRough(Nameclusterfilename, SimScoreRoughNegative, SimScoreRoughPositive)
% check the performance of SimScoreRoughNegative and SimScoreRoughPositive
% to see if the grouping matches with the cluster file.


  Cluster = ReadClusterfile(Nameclusterfilename) ;
  
  %ret = TwofilesInOneCluster(Cluster, file1, file2) ;  
  
  SSLbelief_posi = ReadSimScoreListWithBeliefs(SimScoreRoughPositive) ;
  SSLbelief_nega = ReadSimScoreListWithBeliefs(SimScoreRoughNegative) ;
  
  
  EvaluateSSLwithCluster(SSLbelief_posi, Cluster)  ;
  fprintf('----------------\n');
  EvaluateSSLwithCluster(SSLbelief_nega, Cluster)  ;
  
function EvaluateSSLwithCluster(SSL, Cluster)  

  NumPairs = length( SSL );
  
  for m = 1:NumPairs      
      name1 = SSL(m).name1 ;
      name2 = SSL(m).name2 ;          
      ClusterResult(m) = TwofilesInOneCluster(Cluster, name1, name2) ;      
      SimBelief(m) =   SSL(m).SimBelief ;
      NosimBelief(m) = SSL(m).NosimBelief;      
      
      if 0
      if (ClusterResult(m)==0) & (SimBelief(m)>0.9)
             ShowImg(name1, 1);
             ShowImg(name2, 2);
             
             fprintf('%f %f\n', SimBelief(m), NosimBelief(m));
             pause
      end
      end
      
  end
 
  
  
  for Thr = 0.99:-0.01:0.20
      SimDiff = SimBelief - NosimBelief ;

      Thisres = (SimDiff > Thr) ;

      FalseAccept = sum((ClusterResult==0).*(Thisres)) ;
      FalseReject = sum((ClusterResult).*(Thisres==0)) ;
      CorrectAccept = sum((ClusterResult==1).*(Thisres==1)) ;
      CorrectReject = sum((ClusterResult==0).*(Thisres==0)) ; 
      fprintf('Thr=%f,   FAcce=%d    FReject=%d   CorA=%d    CorRej=%d\n', Thr, FalseAccept, FalseReject, CorrectAccept, CorrectReject);
  end  
   
function ret = TwofilesInOneCluster(Cluster, file1, file2)
% check if two files in the same cluster
  N = length(Cluster);
  pos1 = 0 ;
  pos2 = 0 ;
   for n=1:N       
      K = Cluster(n).size ;        
      for k=1:K          
          filename = Cluster(n).name(k).str ;
          if strcmp(filename, file1)
              pos1 = n ;
              break ;
          end
      end
   end
   for n=1:N       
      K = Cluster(n).size ;        
      for k=1:K          
          filename = Cluster(n).name(k).str ;
          if strcmp(filename, file2)
              pos2 = n ;
              break ;
          end
      end
   end
   
   if pos1==0 | pos2==0
       ret = 0 ;
   else
       if pos1==pos2
           ret = 1 ;
       else
           ret = 0 ;
       end
   end

function ShowImg(filename, idx)
if 0
   I = imread(filename);
   Y = I(:,:,1) ;
   [row, col] = size(Y);
   resizeratio = 512/row ;
   Yt = imresize(Y, resizeratio);   
   figure(idx)
      imshow(uint8(Yt));
      title(filename);
   pos = get(idx, 'position');
   pos(1) = (idx-1)*800 ;
   pos(2) = pos(2) + 100 ;
   set(idx,'position', pos);
else
   figure(idx)
      imshow(imresize(imread(filename), 0.1));
      title(filename);
   pos = get(idx, 'position');
   pos(1) = (idx-1)*750 ;
   pos(2) = pos(2) + 100 ;
   set(idx,'position', pos);
    
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

function SSL = ReadSimScoreListWithBeliefs(filename)
   fin = fopen(filename, 'r');

   if fin==-1
       fprintf('Open sim score list file error\n');
       SSL(1).name1 = 'XXXXX';
       SSL(1).name2 = 'XXXX' ;
       SSL(1).simscore = 0;
       SSL(1).SimBelief = 0 ; % dummy codes
       SSL(1).NosimBelief = 0 ; % dummy codes
       SSL(1).NullBelief = 0 ; % dummy codes  
       %fclose(fin);
       return ;
   end
   
    N = 0 ;
    len = 1;
    while len>0
      [name len] = fscanf(fin, '%s ', 1) ;
      
      if len>0
          
        N = N + 1;
        SSL(N).name1 = name ;
        
        [name len] = fscanf(fin, '%s ', 1) ;
        SSL(N).name2 = name ;
        
        SSL(N).simscore = fscanf(fin, '%f', 1);
        SSL(N).SimBelief = fscanf(fin, '%f', 1);
        SSL(N).NosimBelief = fscanf(fin, '%f', 1);
        SSL(N).NullBelief  = fscanf(fin, '%f', 1);
      end
    end
   
   fclose(fin);
  
  