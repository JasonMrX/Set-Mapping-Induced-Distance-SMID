function compute_allmetrics_forpairs_basedon8(ConsensusPairFile, AllMetrics_SimFile)
% the same as compute_allmetrics_forpairs(), except that the 
% GreedyComputationISmai8 is used instead of GreedyComputationISmai4.
%  Nov. 21. 


% Read pairs info from ConsensusPairFile, compute all metrics and store
% them into AllMetrics_SimFile

%  compute_allmetrics_forpairs('ConsensusPairs.dat', 'AllMetrics_Sim.dat')


 ImagePair = ReadConsensusPairFile(ConsensusPairFile) ;
 
 foutsim=fopen(AllMetrics_SimFile, 'w');

 Npairs = length(ImagePair);
 
 for n=1:Npairs
     % process sim-group
         filename1 = ImagePair(n).name1 ;
         filename2 = ImagePair(n).name2 ;
         
         fprintf('%s %s\n', filename1, filename2);
         
         [success hist_IS mDistPixel mDistDCT pmse] = GreedyComputationISmai8(filename1, filename2) ;
         if success==1 % skip if sizes are not the same 
             [Nmatch1 VScore1  Nmatch2 VScore2] = ComputeSIFT(filename1, filename2, 1.5, 1.5, 0.15) ;            
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
        