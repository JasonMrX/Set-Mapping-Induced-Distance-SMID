function SumVoteResult

  ImageFile = 'TagSys_TestImageInfo.dat';
  VoteFile  = 'TagSys_TestImageVote.dat';
  Namemapfile = 'Int_pairname.dat';

  ImagePair = ReadImageInfo(ImageFile, Namemapfile)  ;
  
  VoteMatrix = ReadVoteMatrix(VoteFile, length(ImagePair));
 
 [Npairs Nvote] = size(VoteMatrix);
 % find pairs with consensus
 ConsensusIndex(1:Npairs) = 0;
 Consensus(1:Npairs) = 0;
 
 Nsim = 0;
 Nnotsim =0;
 Nnotsure =0;
 for n=1:Npairs
   
     meanvote = mean(VoteMatrix(n,:));
     % 1: simimar
     % 2: not-similar
     % 3: not-sure
     ConsensusIndex(n) = 1 ;
     for j=1:Nvote     
         if VoteMatrix(n,j)~=meanvote
             ConsensusIndex(n) = 0 ;
         end
     end
     
     if ConsensusIndex(n)==1
         Consensus(n) = meanvote ;
         if meanvote==1
             Nsim = Nsim+1;
         end
         if meanvote==2
             Nnotsim = Nnotsim+1;
         end
         if meanvote==3
             Nnotsure = Nnotsure+1;
         end
         
     end     
 end
 
 fprintf('Total %d,   Consensused: %d (%d  %d  %d) \n', Npairs, sum(ConsensusIndex), Nsim, Nnotsim, Nnotsure );
 
 
 % write out pairs with consensus
 foutsim=fopen('ConsensusPairs.dat', 'w');

 for imageid=1:Npairs
     % process sim-group
     if ConsensusIndex(imageid)==1 && ( Consensus(imageid)==1 || Consensus(imageid)==2)
         filename1 = ImagePair(imageid).name1 ;
         filename2 = ImagePair(imageid).name2 ;
         
         fprintf(foutsim, '%5d  %16s %16s %16s    %d\n', imageid, ImagePair(imageid).pairname, filename1, filename2,Consensus(imageid));
     end     
 end
 fclose(foutsim);

 %%%%%%%%%%%%%%%%%
 return ;
 


     % illustrate voting 
     VoteMatrix(VoteMatrix==2) = 5 ;
     VoteMatrix(VoteMatrix==3) = 2 ;
     VoteMatrix(VoteMatrix==5) = 3 ;

     for n=1:8
         figure
         plot(1:Npairs, VoteMatrix(:,n) - VoteMatrix(:,9), 'r.-');
         %hold on
         %plot(1:Npairs, VoteMatrix(:,9)+0.05, 'b.-');
         s = sprintf('%d', n);
         title(s);

         %fprintf('%d   %d \n', n, sum(sum(VoteMatrix(:,n)~=VoteMatrix(:,9))));
         %fprintf('     %d \n',  sum(sum(abs(VoteMatrix(:,n)-VoteMatrix(:,9))>1)));
     end

 
  
function VoteMatrix = ReadVoteMatrix(VoteFile, N)
  fin = fopen(VoteFile, 'r');
    
  VoteMatrix = zeros(N, 9); % 8 voters plus one machine results by id=110
  
    
  if fin==-1 
    fprintf('open file error\n');
    return ;
  end
  
  len = 1;
    while len>0
      [userid len] = fscanf(fin, '%d ', 1) ;
      if len>0
       
          if userid>20
              userid = 9 ;
          end
          imageid = fscanf(fin, '%d ', 1);
          vote = fscanf(fin, '%d ', 1);
          
          VoteMatrix(imageid, userid) = vote ;
      end
    end
  
  
  fclose(fin);

function ImagePair = ReadImageInfo(ImageFile, Namemapfile)  
  fin1 = fopen(ImageFile, 'r');
  fin2 = fopen(Namemapfile, 'r');
    
  if fin1==-1  || fin2==-1
    fprintf('open file error\n');
    return ;
  end
  
    
    len = 1;
    while len>0
      [imageid len] = fscanf(fin1, '%d ', 1) ;
      
      if len>0
          len = 1;
    
          pairname = fscanf(fin1, '%s ', 1) ;
          belief = fscanf(fin1, '%f ', 1) ;
          ss = fscanf(fin1, '%f ', 1);
       
          pairnametmp = fscanf(fin2, '%s ', 1);
          name1 = fscanf(fin2, '%s ', 1);
          name2 = fscanf(fin2, '%s ', 1);
          
          if strcmp(pairname, pairnametmp)==0
              fprintf('%s and %s are not synced. please check\n', ImageFile, Namemapfile);
              return ;
          end
          
          ImagePair(imageid).pairname = pairname;
          ImagePair(imageid).name1 = name1;
          ImagePair(imageid).name2 = name2;
          ImagePair(imageid).belief = belief;

          
      end
    end
  
  fclose(fin1);
  fclose(fin2);
  