function ComputeResultDisplay
% display for showing the power of block-matching in pixel domain.

close all;
% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

  [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric('ComputeResult.dat') ;
 
  
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
      xlabel('d_{\phi_{2}^{k, l}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      
      figure,
      low = min([SimMetric(SelectSet1,2); NotsimMetric(SelectSet2,2)]);
      high = max([SimMetric(SelectSet1,2); NotsimMetric(SelectSet2,2)]);
      range2 = low:0.3:high;
      Sim_DCTm = hist(SimMetric(SelectSet1,2), range2);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,2), range2);
      ax = range2;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('Pixel Distance');
      ylabel('Number of pairs of images');
      
      figure,
      SimDist = SimMetric(SelectSet1, 1);
      SimDistCost = SimMetric(SelectSet1, 2);
      NotSimDist = NotsimMetric(SelectSet2,1);
      NotSimDistCost = NotsimMetric(SelectSet2,2);
      plot(SimDist, SimDistCost, 'g.');
      hold on;
      plot(NotSimDist, NotSimDistCost, 'r.');
 
  end
  
  
  
  
function [SimID, NotsimID, SimMetric, NotsimMetric, SimName, NotSimName] = ReadMetric(filename)

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
             SimMetric(pos1, 1:2) = fscanf(fin, '%f ', 2);
             SimID(pos1) = imageid;
             SimName(pos1).name1 = name1 ;
             SimName(pos1).name2 = name2 ;
             
             pos1 = pos1 + 1 ;
         end
         if type ==2
             NotsimMetric(pos2, 1:2) = fscanf(fin, '%f ', 2);
             NotsimID(pos2) = imageid;
             NotSimName(pos2).name1 = name1 ;
             NotSimName(pos2).name2 = name2 ;
             
             pos2 = pos2 + 1 ;
         end
         
          
      end
      
  end
  
  
  fclose(fin);