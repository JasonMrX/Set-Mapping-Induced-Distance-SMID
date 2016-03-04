function ComputeResultDisplay3_ambi
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult3.dat', 5) ;
 
  
  SelectSet1 = (SimID(:)>0)  ;
  SelectSet2 = (NotsimID(:)>0)  ;

  
  if 1
      figure,
      range1 = 0:20:2000*1.25;
      Sim_DCTm    = hist(SimMetric   (SelectSet1,1), range1);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,1), range1);
      ax = range1;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{2, 2}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      simDist = SimMetric(:, 1);
      notSimDist = NotsimMetric(:, 1);
      simDist_ambi = simDist >= min(notSimDist);
      notSimDist_ambi = notSimDist <= max(simDist);
      
%       simName_ambi = SimName(simDist_ambi);
%       simName_ambi(randperm(numel(simName_ambi))) = simName_ambi;
%       for i = 1 : numel(simName_ambi)
%           figure;
%           subplot(1,2,1); imshow(simName_ambi(i).name1);
%           subplot(1,2,2); imshow(simName_ambi(i).name2);
%       end
      
      figure,
      idx = 1;
      low = min([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      high = max([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      range = low: (high - low) / 50: high;
      Sim_DCTm = hist(SimMetric(simDist_ambi,idx), range);
      Notsim_DCTm = hist(NotsimMetric(notSimDist_ambi,idx), range);
      ax = range;
%       plot(ax, Sim_DCTm, 'g-');
%       hold on 
%       plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{2, 2}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      
%       ========================= 8by8 Distortion =========================
%       figure,
%       range1 = 0:20:2000*1.25;
%       Sim_DCTm    = hist(SimMetric   (SelectSet1,5), range1);
%       Notsim_DCTm = hist(NotsimMetric(SelectSet2,5), range1);
%       ax = range1;
%       plot(ax, Sim_DCTm, 'g-');
%       hold on 
%       plot(ax, Notsim_DCTm, 'r-');
%       legend('Similar pairs', 'Not-similar pairs');
%       xlabel('d_{\phi_{2}^{8, 8}} ({\bf A}, {\bf B})');
%       ylabel('Number of pairs of images');
% =============================== P Distance =============================
      
      figure,
      idx = 2;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 100: high;
      Sim_DCTm = hist(SimMetric(SelectSet1, idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2, idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('Pixel Distance');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 2;
      low = min([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      high = max([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      range = low: (high - low) / 50: high;
      Sim_DCTm = hist(SimMetric(simDist_ambi,idx), range);
      Notsim_DCTm = hist(NotsimMetric(notSimDist_ambi,idx), range);
      ax = range;
%       plot(ax, Sim_DCTm, 'g-');
%       hold on 
%       plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('Pixel Distance - ambi');
      ylabel('Number of pairs of images');
%       ================ Variance ======================================
      
      figure,
      idx = 3;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 100: high;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Variance');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 3;
      low = min([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      high = max([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      range = low: (high - low) / 50: high;
      Sim_DCTm = hist(SimMetric(simDist_ambi,idx), range);
      Notsim_DCTm = hist(NotsimMetric(notSimDist_ambi,idx), range);
      ax = range;
%       plot(ax, Sim_DCTm, 'g-');
%       hold on 
%       plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Variance - ambi');
      ylabel('Number of pairs of images');
%      ===================== Entropy ===================================== 
      
      figure,
      idx = 4;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 300: high;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Entropy');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 4;
      low = min([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      high = max([SimMetric(simDist_ambi, idx); NotsimMetric(notSimDist_ambi, idx)]);
      range = low: (high - low) / 50: high;
      Sim_DCTm = hist(SimMetric(simDist_ambi,idx), range);
      Notsim_DCTm = hist(NotsimMetric(notSimDist_ambi,idx), range);
      ax = range;
%       plot(ax, Sim_DCTm, 'g-');
%       hold on 
%       plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Entropy - ambi');
      ylabel('Number of pairs of images');
%       ================ 2D plot =========================================
%       figure,
%       SimDist = SimMetric(SelectSet1, 1);
%       SimDistCost = SimMetric(SelectSet1, 2);
%       NotSimDist = NotsimMetric(SelectSet2,1);
%       NotSimDistCost = NotsimMetric(SelectSet2,2);
%       plot(SimDist, SimDistCost, 'b.');
%       hold on;
%       plot(NotSimDist, NotSimDistCost, 'r.');
%       
%       figure,
%       SimDist = SimMetric(SelectSet1, 1);
%       SimDistCost = SimMetric(SelectSet1, 3);
%       NotSimDist = NotsimMetric(SelectSet2,1);
%       NotSimDistCost = NotsimMetric(SelectSet2,3);
%       plot(SimDist, SimDistCost, 'b.');
%       hold on;
%       plot(NotSimDist, NotSimDistCost, 'r.');
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 4);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,4);
      plot(SimDist, SimDistCost, 'b.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
 
  end
  
  
  
  
function [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric(filename, numMetrics)

  fin = fopen(filename, 'r');
  
  len = 1 ;
  
  pos1 = 1 ;
  pos2 = 1 ;
  while(len>0)
      
      [imageid, len]= fscanf(fin, '%d ', 1);
      
      if len>0
         name1 = fscanf(fin, '%s ', 1);
         name2 = fscanf(fin, '%s ', 1);
         type = fscanf(fin, '%d ', 1);
         
         if type ==1
             SimMetric(pos1, 1:numMetrics) = fscanf(fin, '%f ', numMetrics);
             SimID(pos1) = imageid;
             SimName(pos1).name1 = name1 ;
             SimName(pos1).name2 = name2 ;
             
             pos1 = pos1 + 1 ;
         end
         if type ==2
             NotsimMetric(pos2, 1:numMetrics) = fscanf(fin, '%f ', numMetrics);
             NotsimID(pos2) = imageid;
             NotSimName(pos2).name1 = name1 ;
             NotSimName(pos2).name2 = name2 ;
             
             pos2 = pos2 + 1 ;
         end
         
          
      end
      
  end
  
  
  fclose(fin);