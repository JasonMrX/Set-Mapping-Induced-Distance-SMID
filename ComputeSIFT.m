

function [Nmatch1 VScore1  Nmatch2 VScore2] = ComputeSIFT(file1, file2, match_thr1, match_thr2, DSfactor)

   path(path, 'C:\xiangyu\workspace\SIFT\sift-0.9.19-bin\sift');
   
   [I1 frame1 descrip1] = GetSIFTfeatures(file1, DSfactor);
   [I2 frame2 descrip2] = GetSIFTfeatures(file2, DSfactor);
   
   if 0
       figure
       plotsiftframe(frame1);
       title('frame1');

       figure
       plotsiftframe(frame2);
       title('frame2');

       figure
       plotsiftdescriptor(descrip1);
       title('desc1');

       figure
       plotsiftdescriptor(descrip2);
       title('desc2');
   end
   
   % original matching
   [matches score]= siftmatch(descrip1, descrip2, match_thr1);
   %score   
   if length(score) ==0
       score(1) = 2 ;
   end
   
   fprintf('Num of matches: %d,  mean-score=%f\n', length(score), mean(score));
            if 0
               figure
               plotmatches(I1, I2, frame1, frame2, matches);
               title('match')
            end
   Nmatch1 = length(score);
   VScore1 = mean(score);

   % EMD matching
   [matches2 score2]= siftmatchEMD(descrip1, descrip2, match_thr2);
   fprintf('Num of matches: %d,  mean-score=%f  (EMD_L1)\n', length(score2), mean(score2));
            if 0
               figure
               plotmatches(I1, I2, frame1, frame2, matches2);
               title('match EMD')
            end
   Nmatch2 = length(score2);
   VScore2 = mean(score2);
   

function [I frames descriptors] = GetSIFTfeatures(file, DSfactor)
   I = imread(file);
   
   I = double(rgb2gray(I)/256);
   I = imresize(I, DSfactor);
   [frames, descriptors] = sift(I, 'Verbosity', 0);
   