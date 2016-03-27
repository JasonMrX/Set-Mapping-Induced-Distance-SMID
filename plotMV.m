argc = 5;
[ImagePair, OtherMetric] = ReadDatFile('ComputeResult_moreDat_P2x2_MVE_LB8x8_EMD_HI.dat', argc) ;

for i = 1230 : numel(ImagePair)
    mvfilename = ['mv/', num2str(ImagePair(i).imageid)];
    a = dlmread(mvfilename);
    a1 = a(:, 1 : 2);
    plane1(:, :, 1) = reshape(a1(:, 1), 37, 27)';
    plane1(:, :, 2) = reshape(a1(:, 2), 37, 27)';
    close all;
    figure, 
    subplot(2, 2, 1);
    imshow(plane1(:, :, 1), [-26, 26]);
    subplot(2, 2, 3);
    imshow(plane1(:, :, 2), [-36, 36]);
    subplot(2, 2, 2)
    imshow(ImagePair(i).name1);
    subplot(2, 2, 4);
    imshow(ImagePair(i).name2);
    
    title(num2str(ImagePair(i).imageid));
    
end
    