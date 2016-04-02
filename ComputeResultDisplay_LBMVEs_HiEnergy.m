function ComputeResultDisplay_LBMVEs_HiEnergy
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult_moreDat_LB8x8Org_LBMVEs(HiEnergy).dat', 13) ;
 
%   figure,  
%   hold on;
%   
%   rocCurve = ComputeROC(SimMetric(:, 3), NotsimMetric(:, 3));
%   plot(rocCurve(1, :), rocCurve(2, :), 'g', 'LineWidth', 2); 
%   
%   rocCurve = ComputeROC(SimMetric(:, 4), NotsimMetric(:, 4));
%   plot(rocCurve(1, :), rocCurve(2, :), 'b', 'LineWidth', 2); 
%   
%   rocCurve = ComputeROC(SimMetric(:, 2), NotsimMetric(:, 2));
%   plot(rocCurve(1, :), rocCurve(2, :), 'r', 'LineWidth', 2); 
%   
%   rocCurve = ComputeROC(SimMetric(:, 1), NotsimMetric(:, 1));
%   plot(rocCurve(1, :), rocCurve(2, :), 'k', 'LineWidth', 2); 
% % 
% %   rocCurve = ComputeROC(SimMetric(:, 5), NotsimMetric(:, 5), 2);
% %   plot(rocCurve(1, :), rocCurve(2, :), 'c', 'LineWidth', 2); 
% %   legend('p2x2', 'MVE', 'LB8x8', 'EMD', 'HI');
%   legend('LB8x8', 'EMD', 'MVE', 'p2x2');
%   xlabel('false positive');
%   ylabel('true positive');
%   title('ROC curve');
%   hold off;
  
%   SelectSet1 = (SimID(:)>0)  ;
%   SelectSet2 = (NotsimID(:)>0)  ;
  SelectSet1 = SimMetric(:, 1) > 70;
  SelectSet2 = NotsimMetric(:, 1) < 200;

  
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
      xlabel('LB8x8Org');
      ylabel('Number of pairs of images');
%       ylim([0 50]);
%       xlim([0 2000]);
      for idx = 2 : 13
          figure,
          range = 3 : 0.05 : 6;
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
          rem = mod(idx - 2, 4);
          quo = floor((idx - 2) / 4);
          if rem == 0
              text = 'py';
          elseif rem == 1
              text = 'ey';
          elseif rem == 2
              text = 'px';
          else
              text = 'ex';
          end
          xlabel([num2str(quo), ': ', text]);
          ylabel('Number of pairs of images');
    %       ylim([0 50]);
    %       xlim([0 2000]);
      end
          
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 2);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,2);
      plot(SimDist, SimDistCost, 'g.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
%       
%       figure,
%       SimDist = SimMetric(SelectSet1, 3);
%       SimDistCost = SimMetric(SelectSet1, 2);
%       NotSimDist = NotsimMetric(SelectSet2,3);
%       NotSimDistCost = NotsimMetric(SelectSet2,2);
%       plot(SimDist, SimDistCost, 'g.');
%       hold on;
%       plot(NotSimDist, NotSimDistCost, 'r.');
 
  end
  
  
