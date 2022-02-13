close all;
warning('off','all');
first_len = floor(0.5*len);
Half_len = floor(0.3 *len);

%Function recluster each cluster
fprintf('-----------------------------------------------------------------\n');
fprintf('RECLUSTERING.....\n');   
fprintf('-----------------------------------------------------------------\n');
N_cluster = I(1:first_len); 

for i = 1:length(N_cluster )
    X  = B(i) * (Res/KB); 
    
    if(i < Half_len && X < Max_TADsize)        
      New_cluster_1 = zeros(Limit,1);
      for I = 1:length(index)
          New_cluster_1(index(I,1),1) = cluster(I); % 0.5
      end
       break;
    elseif (X < Max_TADsize)
        
        break;
          
    end
    fprintf('kmeans_Clustering %d ........\n',i);   
    Ind = find(cluster==N_cluster(i)); % Ind contains the rows 
    % find the rows in big data
    Data=  Comb_Data(Ind,:);
    row = size(Data,1);
    Clust = zeros(row);
    % Recluster
    N = size(Data,1);
    No_of_Cluster = floor(sqrt(N/2)); % rough estimatre number of cluster 
    test_num   = 3;    
    [T,~,sumd] = kmeans(Data,No_of_Cluster,'Distance','sqEuclidean','start','uniform', 'emptyaction',...
        'drop');    
    Clust = T;
    destortion_temp=sum(sumd);
    for test_count=2:test_num     
        [V,~,sumd]=kmeans(Data, No_of_Cluster,'Distance','sqEuclidean','emptyaction','drop');
        if ( sum(sumd) < destortion_temp )
           Clust = V;  
           destortion_temp = sum(sumd);
        end
    end 
 
    % Replace thethe rows of the reclustered with new defined cluster
    cluster(Ind) = Clust;
   if(i==Half_len)
           New_cluster_1 = zeros(Limit,1);
          for I = 1:length(index)
               New_cluster_1(index(I,1),1) = cluster(I); % 0.5
          end

    end
end
 % Eventually new clusters are defined 

 
  New_cluster_2 = zeros(Limit,1);
  for I = 1:length(index)
       New_cluster_2(index(I,1),1) = cluster(I); % 0.7
  end 


