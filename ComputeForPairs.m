function ComputeForPairs(AllMetrics_SimFile)
% the same as compute_allmetrics_forpairs(), except that the 
% GreedyComputationISmai8 is used instead of GreedyComputationISmai4.
%  Nov. 21. 
ConsensusPairFile = 'ComputeConsensusPairs.dat';

% Read pairs info from ConsensusPairFile, compute all metrics and store
% them into AllMetrics_SimFile

%  compute_allmetrics_forpairs('ConsensusPairs.dat', 'AllMetrics_Sim.dat')

    ImagePair = ReadConsensusPairFile(ConsensusPairFile) ;

%     Uncomment this part if you want to see a few comparison at the beginning
%     x = randperm(numel(ImagePair));
%     ImagePair(x) = ImagePair;
    
    foutsim=fopen(AllMetrics_SimFile, 'w');

    Npairs = length(ImagePair);

    for n=1:Npairs
     % process sim-group
         filename1 = ImagePair(n).name1 ;
         filename2 = ImagePair(n).name2 ;

         fprintf('%s %s\n', filename1, filename2);

%          [mDistPixel, mDistPixel8, mDistCost, mDistVar, mDistEntropy, success] = ComputeD2Complexity(filename1, filename2)  
         [mDistPixel, mDistEntropy, mDistDCT, success] = ComputeDCTComplexity(filename1, filename2)  

         if success
            fprintf(foutsim, '%d %16s %16s  %d   %7.1f  %7.4f %7.1f\n', ImagePair(n).imageid, filename1, filename2, ImagePair(n).vote, mDistPixel, mDistEntropy, mDistDCT); 
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
        