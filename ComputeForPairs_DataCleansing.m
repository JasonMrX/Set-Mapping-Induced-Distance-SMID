function ComputeForPairs_DataCleansing
% the same as compute_allmetrics_forpairs(), except that the 
% GreedyComputationISmai8 is used instead of GreedyComputationISmai4.
%  Nov. 21. 

% Read pairs info from ConsensusPairFile, compute all metrics and store
% them into AllMetrics_SimFile

%  compute_allmetrics_forpairs('ConsensusPairs.dat', 'AllMetrics_Sim.dat')
    argc = 2;
    [ImagePair, OtherMetric] = ReadMetric('ComputeResult_moreDat_LB8x8Org_LBMVEs(HiEnergy_20).dat', argc) ;
    AllMetrics_SimFile = 'ComputeResult_moreDat_LB8x8Org_LBMVEs(HiEnergy_20)_cleaned.dat';
    vote_valid_file = 'vote_validation.dat';

%     Uncomment this part if you want to see a few comparison at the beginning
%     x = randperm(numel(ImagePair));
%     ImagePair(x) = ImagePair;
%     OtherMetric(x, :) = OtherMetric;
    
    foutsim=fopen(AllMetrics_SimFile, 'a');
    foutvote = fopen(vote_valid_file, 'a');

    Npairs = length(ImagePair);

    format = ['%4d %16s %16s %3d   ', repmat('%10.4f ', 1, 2), '\n'];
    for n = 1405 : Npairs
        filename1 = ImagePair(n).name1 ;
        filename2 = ImagePair(n).name2 ;
        if (ImagePair(n).vote == 1 && OtherMetric(n, 1) <= 70) || (ImagePair(n).vote == 2 && OtherMetric(n, 1) >= 200)
            fprintf(foutsim, format, ImagePair(n).imageid, filename1, filename2, ImagePair(n).vote, OtherMetric(n, :));
            fprintf(foutvote, '%d\n', 1);
            continue;
        end

        fprintf('%d %s %s %d\n', ImagePair(n).imageid, filename1, filename2, ImagePair(n).vote);
        
        sizef = 264*264;
        [ X1oriL] = ReadResizeGrayImage(filename1, sizef);
        [ X2oriL] = ReadResizeGrayImage(filename2, sizef);
        [rowL, colL] = size(X1oriL);
        [row2, col2] = size(X2oriL);
        
        success = 1;
        if rowL ~= row2 || colL ~= col2
            success = 0;
        end 

        if success
            close all;
            figure,
            subplot(1, 2, 1);
            imshow(uint8(X1oriL));
            title(OtherMetric(n, 1));
            subplot(1, 2, 2);
            imshow(uint8(X2oriL));
            title(OtherMetric(n, 2));
            valid = input('valid? ');
            fprintf(foutvote, '%d\n', valid);
            if valid
                fprintf(foutsim, format, ImagePair(n).imageid, filename1, filename2, ImagePair(n).vote, OtherMetric(n, :));
            end
        end

    end
    fclose(foutsim);
    fclose(foutvote);

    return ;
        
function [ImagePair, Metric] = ReadMetric(filename, argc)
  fin = fopen(filename, 'r');
  len = 1 ;
  pos1 = 1 ;
  pos2 = 1 ;
  while (len > 0)
      [imageid len] = fscanf(fin, '%d ', 1);
      if len > 0
         name1 = fscanf(fin, '%s ', 1);
         name2 = fscanf(fin, '%s ', 1);
         type = fscanf(fin, '%d ', 1);
         ImagePair(pos1).imageid = imageid;
         
         ImagePair(pos1).name1 = name1;
         ImagePair(pos1).name2 = name2;
         ImagePair(pos1).vote = type;
         Metric(pos1, 1 : argc) = fscanf(fin, '%f ', argc);
         pos1 = pos1 + 1 ;
      end
  end
  fclose(fin);
  