function DisplayAllMetric_basedon10
% display for showing the power of block-matching in pixel domain.

% use 9 as base. 10 uses 2x2 DCT blocksize

% 1  Histogram
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2
  [SimID  NotsimID  SimMetric NotsimMetric  SimName NotSimName] = ReadMetric('testresult10.dat') ;
  
  SelectSet1 = (SimID(:)>0)  ;
  SelectSet2 = (NotsimID(:)>0)  ;
  
  if 1
      figure
      Sim_DCTm = hist(SimMetric(SelectSet1,3), 0:5:2000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,3), 0:5:2000);
      ax = 0:5:2000;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('T_d(Q_A^{2,2}, Q_B^{2,2})');
      ylabel('Number of pairs of images');
      
      figure
      Sim_DCTm    = hist(SimMetric   (SelectSet1,2), 0:40:3000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,2), 0:40:3000);
      plot(Sim_DCTm, 'g-');
      hold on 
      plot(Notsim_DCTm, 'r--');

      figure
      Sim_DCTm = hist(SimMetric(SelectSet1,1), 0:0.01:1);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,1), 0:0.01:1);
      plot(Sim_DCTm, 'g-');
      hold on 
      plot(Notsim_DCTm, 'r--');
    else
      
      figure
      ax =  0:5:2500 ;
      Sim_DCTm = hist(SimMetric(SelectSet1,3), 0:5:2500);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,3), 0:5:2500);
      bar(ax, Sim_DCTm, 'b');
      hold on 
      bar(ax, Notsim_DCTm, 'r');


      figure
      Sim_DCTm    = hist(SimMetric   (SelectSet1,2), 0:40:3000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,2), 0:40:3000);
      ax = 0:40:3000 ;
      bar(ax, Sim_DCTm, 'b');
      hold on 
      bar(ax, Notsim_DCTm, 'r');

      figure
      Sim_DCTm = hist(SimMetric(SelectSet1,1), 0:0.01:1);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,1), 0:0.01:1);
      ax = 0:0.01:1 ;
      bar(ax, Sim_DCTm, 'b');
      hold on 
      bar(ax, Notsim_DCTm, 'r');
    
  end
  
  
  Sim_MatchN1 = SimMetric(SelectSet1, 5) ;
  Sim_Score1 = SimMetric(SelectSet1, 6) ;
  Sim_MatchN2 = SimMetric(SelectSet1, 7) ;
  Sim_Score2 = SimMetric(SelectSet1, 8) ;
  Notsim_MatchN1 = NotsimMetric(SelectSet2, 5) ;
  Notsim_Score1  = NotsimMetric(SelectSet2, 6) ;
  Notsim_MatchN2 = NotsimMetric(SelectSet2, 7) ;
  Notsim_Score2  = NotsimMetric(SelectSet2, 8) ;
  figure
  plot(Sim_MatchN1, Sim_Score1, 'b.');
  hold on 
  plot(Notsim_MatchN1, Notsim_Score1, 'rd');
  
  figure
  plot(Sim_MatchN2, Sim_Score2, 'b.');
  hold on 
  plot(Notsim_MatchN2, Notsim_Score2, 'rd');
  
  fprintf('NumSim=%d   Num_Notsim=%d \n', sum(SelectSet1), sum(SelectSet2));
  
  
  
  
  
function [SimID  NotsimID  SimMetric NotsimMetric  SimName NotSimName] = ReadMetric(filename)

  fin = fopen(filename, 'r');
  
  len = 1 ;
  
  pos1 = 1 ;
  pos2 = 1 ;
  while(len>0)
      
      [imageid len]= fscanf(fin, '%d ', 1);
      
      if len>0
         name1 = fscanf(fin, '%s ', 1);
         name2 = fscanf(fin, '%s ', 1);
         type = fscanf(fin, '%d ', 1);
         
         if type ==1
             SimMetric(pos1, 1:8) = fscanf(fin, '%f ', 8);
             SimID(pos1) = imageid;
             SimName(pos1).name1 = name1 ;
             SimName(pos1).name2 = name2 ;
             
             pos1 = pos1 + 1 ;
         end
         if type ==2
             NotsimMetric(pos2, 1:8) = fscanf(fin, '%f ', 8);
             NotsimID(pos2) = imageid;
             NotSimName(pos2).name1 = name1 ;
             NotSimName(pos2).name2 = name2 ;
             
             pos2 = pos2 + 1 ;
         end
         
          
      end
      
  end
  
  
  fclose(fin);