
warning('off','all');
first_len = floor(0.5*len);
Half_len = floor(0.3 *len);
explore_limit = 3;
%Function recluster each cluster
fprintf('-----------------------------------------------------------------\n');
fprintf('RECLUSTERING.....\n');   
fprintf('-----------------------------------------------------------------\n');
N_cluster = I(1:first_len); 

for i = 1:length(N_cluster )
    X  = B(i) * (Res/KB); % minimum is approximately 200kb
    if (X < 600)
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
    k_opt =  No_of_Cluster;
    Start = k_opt - explore_limit;
    End = k_opt + explore_limit;
    if (Start < 2)
        Start =    k_opt;
    end    
    leng = End - Start + 1;
    test_num   = 2; 
    local_Clusters = zeros(row,leng);
    K = Start;
    for j = 1:leng
        No_of_Cluster = K;   % No of clusters        
        [T,~,sumd] = kmeans(Data,No_of_Cluster,'Distance','sqEuclidean','start','uniform', 'emptyaction',...
            'drop');    
        local_Clusters(:,j) = T;
        destortion_temp=sum(sumd);
        for test_count=1:test_num     
            [V,~,sumd]=kmeans(Data, No_of_Cluster,'Distance','sqEuclidean','emptyaction','drop');
            if ( sum(sumd) < destortion_temp )
              local_Clusters(:,j) = V;  
               destortion_temp = sum(sumd);
            end
        end
        K = K + 1;
        Agn = Order_TADNum( local_Clusters(:,j)); % Order the cluster
        local_Clusters(:,j) =  Agn';
    end
    %---------------------------------------------------------------------
    % Extract TD from Cluster and check Quality
    %---------------------------------------------------------------------
    upper_limit = Ind(1) + length(Ind) -1 ;
    Initial_diff = -50;  
    Diff=   [];
    Base_index = Ind(1) -1 ;
    for j = 1:leng        
        [TD] = F_TAD_recluster(Chr_Data, local_Clusters(:,j),Base_index); 
        QualityCheck_recluster;
        Diff = [Diff; Intra_Inter_diff];
        if ( Intra_Inter_diff >=  Initial_diff )
            ind = j; % Index with highest difference
            Initial_diff = Intra_Inter_diff;
        end        
    end
    fprintf('chosen = %d\n ',ind);
    close all;
    %---------------------------------------------------------------------
 
    % Replace thethe rows of the reclustered with new defined cluster  
    cluster(Ind) =  local_Clusters(:,ind);
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



