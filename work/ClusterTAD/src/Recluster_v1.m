%% Find TAD from the result

C = dlmread(clustfile);
Clusters = [];
col  =length(C(1,:));
for k = 1:col
    c = C(:,k);
    cluster = [];
    for i = 1:length(C)    
        T = find(zeroRows==i);
        if (isempty(T)==1)
            cluster = [cluster; c(i)]; 
        end
    end

    % Resiign cluster label
     Agn = Order_TADNum(cluster);
     cluster =  Agn';


    % Get the size of each cluster
    S = [];
    for i = 1:max(cluster)
    t = find(cluster(:,1)==i);	
        S = [S;length(t)];
    end
    % Rank the cluster based on the size
    [B,I] = sort(S,'descend');
    % Save Cluster assignment Alone
    foldname = [Resultpath,'/Re-Clusters'];
    if ~exist(foldname, 'dir')
        % Folder does not exist so create it.
        mkdir(foldname);
    end
    out_path_1 = [foldname,'/'];
    clustfile_new=strcat(out_path_1,algorithm,'_Cluster_',name,'_.txt');

    len  = length(B);
    % pick the First half OR LESS
    %-------------------------------------------------------------------------
     reclust_v1; % RECLUSTER
    try        
         Clusters = [Clusters New_cluster_1,New_cluster_2 ];    
    catch ME
        disp('===================== RECLUSTERING  Error Occured ==========================');
        disp('Reclustering  NOT REQUIRED for this Dataset ');
        disp('Find the name of the best TAD in the Readmefile in Quality folder ../Quality/');
        disp( msg);
        disp('Find the TADs identified in the TAD/ folder');
         disp('============================================================================');
        break;       
    end
    %------------------------------------------------------------------------    
end

dlmwrite(clustfile_new,Clusters);
Quality_Assessment_v2; %Quality Assessment
QualityCheck_v2; %Quality Check