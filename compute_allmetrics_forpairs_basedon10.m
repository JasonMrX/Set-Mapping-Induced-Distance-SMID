function compute_allmetrics_forpairs_basedon10(ConsensusPairFile, AllMetrics_SimFile)
% the same as compute_allmetrics_forpairs_basedon9(), except that the 
% GreedyComputationISmai10 is used instead of GreedyComputationISmai9.
%  Dec. 10  

% Note GreedyComputationISmai10() only calculate DCT-based metric.
% Other metric are read from AllMetrics_Sim_for9.dat directly.

% compute_allmetrics_forpairs_basedon10('ConsensusPairs.dat', 'AllMetrics_Sim_for10.dat')
 
% compare with AllMetrics_Sim_for9.dat, AllMetrics_Sim_for10.dat only
% changes the DCT-based metric to use 2x2 block size.

 %ImagePair = ReadConsensusPairFile(ConsensusPairFile) ;
 
 foutsim=fopen(AllMetrics_SimFile, 'w');

 [ImagePair, OtherMetric] = ReadMetric('testresult9.dat') ;
 
 Npairs = length(ImagePair);
 
 pos = 1; 
 
 for n=1:Npairs
     % process sim-group
         filename1 = ImagePair(n).name1 ;
         filename2 = ImagePair(n).name2 ;
     
          pos = n;
         [success  mDistDCT ] = GreedyComputationISmai10(filename1, filename2) ;
         if success==1 % skip if sizes are not the same 
                      fprintf('%s %s  T8(%8.3f)  T2(%8.3f)\n', filename1, filename2, OtherMetric(pos, 3), mDistDCT );
                      
             hist_IS = OtherMetric(pos, 1);
             mDistPixel = OtherMetric(pos,2); 
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
  
        