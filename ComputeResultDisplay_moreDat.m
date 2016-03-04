function ComputeResultDisplay_moreDat
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult_moreDat.dat', 5) ;
 
  
  SelectSet1 = (SimID(:)>0)  ;
  SelectSet2 = (NotsimID(:)>0)  ;

  
  if 1
      figure,
      range1 = 0:20:2000*1.25;
      Sim_DCTm    = hist(SimMetric   (SelectSet1,1), range1);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,1), range1);
      ax = range1;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{2, 2}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      
      figure,
      range1 = 0:20:2000*1.25;
      Sim_DCTm    = hist(SimMetric   (SelectSet1,5), range1);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,5), range1);
      ax = range1;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2}^{8, 8}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');

      figure,
      idx = 2;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 100: high;
      Sim_DCTm = hist(SimMetric(SelectSet1, idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2, idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('Pixel Distance');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 3;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 100: high;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Variance');
      ylabel('Number of pairs of images');
      
      figure,
      idx = 4;
      low = min([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      high = max([SimMetric(SelectSet1, idx); NotsimMetric(SelectSet2, idx)]);
      range = low: (high - low) / 300: high;
      Sim_DCTm = hist(SimMetric(SelectSet1,idx), range);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,idx), range);
      ax = range;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('MV Entropy');
      ylabel('Number of pairs of images');
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 2);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,2);
      plot(SimDist, SimDistCost, 'b.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 3);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,3);
      plot(SimDist, SimDistCost, 'b.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 4);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,4);
      plot(SimDist, SimDistCost, 'g.');
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