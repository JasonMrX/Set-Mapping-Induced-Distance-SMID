function ShowSomePairs

[SimID  NotsimID  SimMetric NotsimMetric  SimName NotSimName] = ReadMetric('ShowSomePairs.dat') ;

SimMetrictmp = SimMetric ;
[Y I] = sort( SimMetrictmp(:,2) );

SimMetric = SimMetric(I, :);
SimName = SimName(I);

      [row col]= size(SimMetric);
      for n=1:row
        fprintf('%s %s    %f  %f  %f  %f  %f  %f  %f  %f \n', SimName(n).name1, SimName(n).name2, SimMetric(n,1), SimMetric(n,2), SimMetric(n,3), SimMetric(n,4), SimMetric(n,1+4), SimMetric(n,2+4), SimMetric(n,3+4), SimMetric(n,4+4));
        s1=sprintf('%d: %5.3f [%5.1f] (%5.1f) %5.1f \n', n,  SimMetric(n,1), SimMetric(n,2), SimMetric(n,3), SimMetric(n,4));
        s2=sprintf('(%4.1f %4.3f) (%4.1f %4.3f)\n', SimMetric(n,1+4), SimMetric(n,2+4), SimMetric(n,3+4), SimMetric(n,4+4));
        figure
        subplot(1,2,1);
        bwimshow(SimName(n).name1);
        title(s1);
        subplot(1,2,2);
        bwimshow(SimName(n).name2);    
        title(s2);
      end

      [row col]= size(NotsimMetric);
      for n=1:row
        fprintf('%s %s    %f  %f  %f  %f  %f  %f  %f  %f \n', NotSimName(n).name1, NotSimName(n).name2, NotsimMetric(n,1), NotsimMetric(n,2), NotsimMetric(n,3), NotsimMetric(n,4), NotsimMetric(n,1+4), NotsimMetric(n,2+4), NotsimMetric(n,3+4), NotsimMetric(n,4+4));
        s=sprintf('Not %s %s %5.3f  %5.1f  %5.1f  %5.1f\n', NotSimName(n).name1, NotSimName(n).name2, NotsimMetric(n,1), NotsimMetric(n,2), NotsimMetric(n,3), NotsimMetric(n,4));
        figure
        subplot(1,2,1);
        bwimshow(NotSimName(n).name1);
        subplot(1,2,2);
        bwimshow(NotSimName(n).name2);    
        title(s);
      end


function bwimshow(filename) 
  I = imread(filename);   
  Y= rgb2gray(I);
  imshow(Y);
  
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
         
         if type ~=1 & type ~=2
             x = fscanf(fin, '%f ', 8); % throw away 
         end
         
         
         
          
      end
      
  end
  
  
  fclose(fin);
  