function ComputeResultDisplay_moreDat_DCT
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult_moreDat_DCT.dat', 3) ;
 
  
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
      c = stem(ax, [Sim_DCTm', Notsim_DCTm']);
      set(c(1), 'Color', 'b');
      set(c(2), 'Color', 'r');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{2, 2}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      ylim([0 50]);
      xlim([0 2000]);

      figure,
      idx = 2;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 200: high;
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
      range = low: (high - low) / 400: high;
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
      ylim([0 50]);
      xlim([0 2000]);
      
   
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
  
  
