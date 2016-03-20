function ComputeResultDisplay_moreDat_P2x2_MVE_LB8x8_EMD_HI
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult_moreDat_P2x2_MVE_LB8x8_EMD_HI.dat', 5) ;
 
  figure,  
  hold on;
  
  rocCurve = ComputeROC(SimMetric(:, 3), NotsimMetric(:, 3));
  plot(rocCurve(1, :), rocCurve(2, :), 'g', 'LineWidth', 2); 
  
  rocCurve = ComputeROC(SimMetric(:, 4), NotsimMetric(:, 4));
  plot(rocCurve(1, :), rocCurve(2, :), 'b', 'LineWidth', 2); 
  
  rocCurve = ComputeROC(SimMetric(:, 2), NotsimMetric(:, 2));
  plot(rocCurve(1, :), rocCurve(2, :), 'r', 'LineWidth', 2); 
  
  rocCurve = ComputeROC(SimMetric(:, 1), NotsimMetric(:, 1));
  plot(rocCurve(1, :), rocCurve(2, :), 'k', 'LineWidth', 2); 
% 
%   rocCurve = ComputeROC(SimMetric(:, 5), NotsimMetric(:, 5), 2);
%   plot(rocCurve(1, :), rocCurve(2, :), 'c', 'LineWidth', 2); 
%   legend('p2x2', 'MVE', 'LB8x8', 'EMD', 'HI');
  legend('LB8x8', 'EMD', 'MVE', 'p2x2');
  xlabel('false positive');
  ylabel('true positive');
  title('ROC curve');
  hold off;
  
  SelectSet1 = (SimID(:)>0)  ;
  SelectSet2 = (NotsimID(:)>0)  ;

  
  if 1
      figure,
      idx = 1;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 400: high;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
%       c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
%       set(c(1), 'Color', 'b');
%       set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{2, 2}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
%       ylim([0 50]);
%       xlim([0 2000]);
        

      figure,
      idx = 2;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 400: high;
      Sim_DCTm = hist(SimMetric(SelectSet1, idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2, idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Entropy');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 3;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
%       range = low: (high - low) / 400: high;  
      range = low : 20 : high;
%       range = 0:5:2000;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('Lower Bound');
      ylabel('Number of pairs of images');
      ylim([0 220]);
      xlim([0 2400]);
      
      figure,
      idx = 4;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
%       range = low: (high - low) / 400: high;
      range = low : 20 : high; 
%       range = 0:5:2000;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('EMD');
      ylabel('Number of pairs of images');
      ylim([0 220]);
      xlim([0 2400]);
      
      figure,
      idx = 5;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 100: high;
%       range = 0:5:2000;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'b-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
%       c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
%       set(c(1), 'Color', 'b');
%       set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('HI');
      ylabel('Number of pairs of images');
%       ylim([0 50]);
%       xlim([0 2000]);
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 2);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,2);
      plot(SimDist, SimDistCost, 'g.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
      
      figure,
      SimDist = SimMetric(SelectSet1, 3);
      SimDistCost = SimMetric(SelectSet1, 2);
      NotSimDist = NotsimMetric(SelectSet2,3);
      NotSimDistCost = NotsimMetric(SelectSet2,2);
      plot(SimDist, SimDistCost, 'g.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
 
  end
  
  
