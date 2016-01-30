function Generate_SimScoreFile_reference(clusterfile, SimScoreRoughPositiveFile, SimScoreRoughNegativeFile, SSLoutfile)
% generate a SSL file, (a version to be viewed) to be used for reference for training and testing




% generate three sets and then get their union
    % Set1: SimPair in RoughPositive file    
    % Set2: all pairs within each cluster, that dont appear in Set1
    % Set3: Selected SimPair in RoughNegative file
    SSL_all_positive = ReadSimScoreListWithBeliefs(SimScoreRoughPositiveFile) ;
    
    SSL_pairs = FromClusterToInnerPair(clusterfile, SSL_all_positive) ;% all inner pairs except that appears in the allpositive set
    %inner pairs here dont have simscores yet. 
    
    SSL_all_negative = ReadSimScoreListWithBeliefs(SimScoreRoughNegativeFile) ;
    [SSL_selected_negative SSL_pairs] = SelectNegativePairs(SSL_all_negative, SSL_pairs);
    
    fout = fopen(SSLoutfile, 'w');
    writeOneSetToFile(fout, SSL_all_positive) ; % all positive pairs
    writeOneSetToFile(fout, SSL_pairs) ;        % other inners pairs, some are initiated from the negative set
    writeOneSetToFile(fout, SSL_selected_negative) ; % FalseAccept in the negative set.   
    fclose(fout);

    fprintf('size of SSL_positive: %d \n', length(SSL_all_positive));
    fprintf('size of other inner pairs: %d \n', length(SSL_pairs));
    fprintf('size of selected negative pairs with false accept: %d \n', length(SSL_selected_negative));
    
function writeOneSetToFile(fout, SSL)    
  Npair = length(SSL);
  
  for m=1:Npair
    fprintf(fout, '%15s  %15s  %5.3f  %5.3f  %5.3f  %5.3f\n', SSL(m).name1, SSL(m).name2, SSL(m).simscore, ...
        SSL(m).SimBelief, SSL(m).NosimBelief, SSL(m).NullBelief);
    
  end
    %Note: Negative set has applied 
    % beliefs[0]>0.08 || beliefs[2]>beliefs[1], as in cluster.cpp
function [SSL_selected_negative  SSL_pairsOut] = SelectNegativePairs(SSL_all_negative, SSL_pairsIn)
% SSL_pairsIn: all inner pairs except that in SimScorePositive; no scores set yet.
% SSL_pairsOut:  scores are set if found in the SimScoreNegative file

% SSL_selelcted_negative: contains some falseAccept pairs in the negative
% set, under a condition such as SimBelief - NosimBelief>0.7
  Npair = length(SSL_all_negative);
  Nselect = 0 ;
  
  SSL_pairsOut = SSL_pairsIn ;
  
  for m=1:Npair
  
      if SSL_all_negative(m).SimBelief - SSL_all_negative(m).NosimBelief>0.7  
          % correctAccept + FalseAccept
          pos = CheckApairInSSL(SSL_pairsIn, SSL_all_negative(m).name1, SSL_all_negative(m).name2);
          if(pos>0)  % CorrectAccept
              SSL_pairsOut(pos) = SSL_all_negative(m);
          else       % falseAccept  Conditonally selected
              Nselect = Nselect + 1 ;
              SSL_selected_negative(Nselect) = SSL_all_negative(m); 
          end
      else
          % get FalseReject
          pos = CheckApairInSSL(SSL_pairsIn, SSL_all_negative(m).name1, SSL_all_negative(m).name2);
          if(pos>0)  % get the reference beliefs 
              SSL_pairsOut(pos) = SSL_all_negative(m);
          end
          % CorrectReject are too many to be included.
      end
  end


% form inner paris for each cluster, but exclude whatever appears in
% SSL_positive
function SSL_pairs = FromClusterToInnerPair(clusterfile, SSL_positive)
  Cluster = ReadClusterfile(clusterfile);
  pairN = 0 ;
  N = length(Cluster);   
   for n=1:N       
      K = Cluster(n).size ;                 
      for k=1:K                    
          filename_1 = Cluster(n).name(k).str ;
          
          for j=k+1:K
              filename_2 = Cluster(n).name(j).str ;
              
              pos_ref = CheckApairInSSL(SSL_positive, filename_1, filename_2) ;
                  
              if (pos_ref>0)                  
                  % fprintf('RefScore is (%5.3f  %5.3f  %5.3f)\n', SSL_ref(pos_ref).SimBelief,SSL_ref(pos_ref).NosimBelief,SSL_ref(pos_ref).NullBelief );
              else
                  pairN = pairN + 1 ;
                  SSL_pairs( pairN ).name1 = filename_1 ;
                  SSL_pairs( pairN ).name2 = filename_2 ;
                  SSL_pairs( pairN ).simscore = -1 ; % 
                  SSL_pairs( pairN ).SimBelief = -1 ; % 
                  SSL_pairs( pairN ).NosimBelief = -1 ; % 
                  SSL_pairs( pairN ).NullBelief = -1 ; % 
              end                  
          end              
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

% pos ==0: not found
% pos is the position in SSL
function pos = CheckApairInSSL(SSL, filename1, filename2)
  
  NumPairs = length( SSL );
  pos = 0 ;
  
  for m = 1:NumPairs
      
      name1 = SSL(m).name1 ;
      name2 = SSL(m).name2 ;    
      
      if strcmp(name1, filename1) 
          if strcmp(name2, filename2)
              pos = m;
              return ;
          end
      end
      if strcmp(name1, filename2) 
          if strcmp(name2, filename1)
              pos = m ;
              return ;
          end
      end      
  end

     
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
      
      %if N>5000
      %    break;
      %end
    end
   
   fclose(fin);
  
