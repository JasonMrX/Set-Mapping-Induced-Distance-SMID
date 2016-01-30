function DisplayAllMetric_basedon13

% copy from DisplayAllMetric_basedon8, but modified to show figures
% for showing the first/second terms as 
%     1st       hist_IS = hist_EMD;
%     2nd        mDistPixel = hist_EMD_shiftinv; 

% 

% 1  Autocorrologram distance
% 2: mDistPixel, 
% 3: mDistDCT, 
% 4: pmse, Nmatch1, VScore1,  Nmatch2, VScore2

[SimID9  NotsimID9  SimMetric9 NotsimMetric9  SimName9 NotSimName9] = ReadMetric('AllMetrics_Sim_for9.dat') ;
% The 9 file is read to do scanning

  [SimID  NotsimID  SimMetric NotsimMetric  SimName NotSimName] = ReadMetric('AllMetrics_Sim_for13.dat') ;
  
 
  
  SelectSet1 = (SimMetric9(:,3)<80)      & (SimMetric9(:,2)<400)  ;
  SelectSet2 = (NotsimMetric9(:,3)>120)  & (NotsimMetric9(:,2)>440)  ;
  
  if 1
      
      figure
      Sim_DCTm = hist(SimMetric(SelectSet1,3), 0:5:2000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,3), 0:5:2000);
      ax =  0:5:2000 ;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('T_d(Q_A^{k,l}, Q_B^{k,l})');
      ylabel('Number of pairs of images');


      figure
      Sim_DCTm    = hist(SimMetric   (SelectSet1,1), 0:24:2000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,1), 0:24:2000);
      ax = 0:24:2000 ;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2, \tau}^{4, 4}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');

      figure
      Sim_DCTm    = hist(SimMetric   (SelectSet1,2), 0:24:2000);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,2), 0:24:2000);
      ax = 0:24:2000 ;
      plot(ax, Sim_DCTm, 'g-');
      hold on 
      plot(ax, Notsim_DCTm, 'r-');
      legend('Similar pairs', 'Not-similar pairs');
      xlabel('d_{\phi_{2, \tau}^{8, 8}} ({\bf A}, {\bf B})');
      ylabel('Number of pairs of images');
      
 
      
    else
      
      figure
      Sim_DCTm = hist(SimMetric(SelectSet1,3), 0:5:2500);
      Notsim_DCTm = hist(NotsimMetric(SelectSet2,3), 0:5:2500);
      ax =  0:5:2500 ;
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