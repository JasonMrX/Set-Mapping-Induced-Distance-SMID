function compute_allmetrics_forpairs_basedon13(ConsensusPairFile, AllMetrics_SimFile)
% the same as compute_allmetrics_forpairs_basedon12(), except that the 
% GreedyComputationISmai13 is used instead of GreedyComputationISmai12.
% To compute blockmatching distance using 4x4 blocksize and 8x8 blocksize
%  Dec. 20, for computing autocorrelogram metric

% Note GreedyComputationISmai13() only calculate 4x4 and 8x8 block matching
% distance.
% Other metric are read from AllMetrics_Sim_for9.dat directly.

% use:
% compute_allmetrics_forpairs_basedon13('ConsensusPairs.dat', 'AllMetrics_Sim_for13.dat')

 
 foutsim=fopen(AllMetrics_SimFile, 'w');

 [ImagePair, OtherMetric] = ReadMetric('testresult9.dat') ;
 
 Npairs = length(ImagePair);
 
 pos = 1; 
 
 for n=1:Npairs
     % process sim-group
         filename1 = ImagePair(n).name1 ;
         filename2 = ImagePair(n).name2 ;
     
          pos = n;
         [success   d4x4 d8x8] = GreedyComputationISmai13(filename1, filename2) ;
         if success==1 % skip if sizes are not the same 
                      %fprintf('%s %s  T8(%8.3f)  T2(%8.3f)\n', filename1, filename2, OtherMetric(pos, 3), mDistDCT );
                      fprintf('%s %s \n', filename1,filename2);
             hist_IS = d4x4;
             mDistPixel = d8x8; 
             mDistDCT = OtherMetric(pos,3) ;
             pmse = OtherMetric(pos,4) ;
             Nmatch1 = OtherMetric(pos,5) ;
             VScore1 = OtherMetric(pos,6) ;
             Nmatch2 = OtherMetric(pos,7) ;
             VScore2 = OtherMetric(pos,8) ;
             %pos = pos + 1 ;
             %[Nmatch1 VScore1  Nmatch2 VScore2] = ComputeSIFT(filename1, filename2, 1.5, 1.5, 0.15) ;            
             fprintf(foutsim, '%d %16s %16s  %d   %6.5f  %7.1f  %7.1f  %7.1f    %d %5.3f   %d %5.3f\n', ImagePair(n).imageid, filename1, filename2, ImagePair(n).vote, hist_IS, mDistPixel, mDistDCT, pmse, Nmatch1, VScore1,  Nmatch2, VScore2); 
         end
     
 end
 fclose(foutsim);

 return ;
 
function ImagePair = ReadConsensusPairFile(filename)
        fin = fopen(filename);
        
        len=1;
        pos = 1 ;
        while len>0
           [imageid len] = fscanf(fin, '%d ', 1);
           if len>0
              ImagePair(pos).imageid = imageid;
              ImagePair(pos).pairname = fscanf(fin, '%s ', 1);
              ImagePair(pos).name1 = fscanf(fin, '%s ', 1);
              ImagePair(pos).name2 = fscanf(fin, '%s ', 1);
              ImagePair(pos).vote = fscanf(fin, '%d' , 1); 
              pos = pos + 1;
           end
            
        end
        fclose(fin);
        
function [ImagePair, Metric] = ReadMetric(filename)

  fin = fopen(filename, 'r');
  
  len = 1 ;
  
  pos1 = 1 ;
  pos2 = 1 ;
  while(len>0)
      
      [imageid len]= fscanf(fin, '%d ', 1);
      
      if len>0
         name1 = fscanf(fin, '%s ', 1);
         name2 = fscanf(fin, '%s ', 1);
         type = fscanf(fin, '%d ', 1);
         ImagePair(pos1).imageid = imageid;
         
         ImagePair(pos1).name1 = name1;
         ImagePair(pos1).name2 = name2;
         ImagePair(pos1).vote = type;
             Metric(pos1, 1:8) = fscanf(fin, '%f ', 8);
             pos1 = pos1 + 1 ;
         
          
      end
      
  end
  
  
  fclose(fin);
  
        