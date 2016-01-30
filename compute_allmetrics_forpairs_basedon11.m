function compute_allmetrics_forpairs_basedon11(ConsensusPairFile, AllMetrics_SimFile)
% the same as compute_allmetrics_forpairs_basedon10(), except that the 
% GreedyComputationISmai11 is used instead of GreedyComputationISmai10.
%  Dec. 12 
 
% Note GreedyComputationISmai11() only calculate hist-EMD with shift
% invariance
% Other metric are read from AllMetrics_Sim_for9.dat directly.

% use:
% compute_allmetrics_forpairs_basedon11('ConsensusPairs.dat', 'AllMetrics_Sim_for11.dat')

 
 foutsim=fopen(AllMetrics_SimFile, 'w');

 [ImagePair, OtherMetric] = ReadMetric('testresult9.dat') ;
 
 Npairs = length(ImagePair);
 
 pos = 1; 
 
 for n=1:Npairs
     % process sim-group
         filename1 = ImagePair(n).name1 ;
         filename2 = ImagePair(n).name2 ;
     
          pos = n;
         [success   hist_EMD hist_EMD_shiftinv ] = GreedyComputationISmai11(filename1, filename2) ;
         if success==1 % skip if sizes are not the same 
                      %fprintf('%s %s  T8(%8.3f)  T2(%8.3f)\n', filename1, filename2, OtherMetric(pos, 3), mDistDCT );
                      fprintf('%s %s \n', filename1,filename2);
             hist_IS = hist_EMD;
             mDistPixel = hist_EMD_shiftinv; 
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
  
        