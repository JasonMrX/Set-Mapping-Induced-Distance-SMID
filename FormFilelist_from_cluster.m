function FormFilelist_from_cluster(ClusterFileIn, namelistfileOut)
% FormFilelist_from_cluster('test_cluster.dat', 'test_groupedfile_list.dat')

   Cluster = ReadClusterfile(ClusterFileIn) ;
    
   N = length(Cluster);   
   filenum = 0 ;   
   for n=1:N       
      K = Cluster(n).size ;                 
      for k=1:K 
          
          thisname = Cluster(n).name(k).str ;
          if filenum>0  % check for duplications
              for z=1:filenum
                  if strcmp( filelist(z), thisname ) 
                      fprintf('Some file appears in more than one cluster. Wrong\n');
                      return ;
                  end
              end
          end          
          
          filenum = filenum + 1 ;                    
          filelist( filenum).name = thisname ;
      end
   end
   
   % write out a filelist 
   fout = fopen(namelistfileOut, 'w');
   fprintf(fout, '%d\n', filenum);
   for n=1:filenum
       fprintf(fout, '%s\n', filelist(n).name);
   end
   fclose(fout);





function Cluster = ReadClusterfile(filename)
  fin = fopen(filename);
  
   if ~fin
       fprintf('Creating a new cluster file\n');
       
       fclose(fin);
       return ;
       
   end
   
    K = 0 ;
    count = 1;
    while count>0
      [M count] = fscanf(fin, '%d', 1) ;
      
      if count>0
          
        K = K + 1;
        
        Cluster(K).size = M ;        
        for m=1:M
          [name len] = fscanf(fin, '%s', 1) ;
          if len>0
              Cluster(K).name(m).str = name ;
          else
              fprintf('reading cluster file error\n');
          end
        end
        
      end
    end
  fclose(fin);      
