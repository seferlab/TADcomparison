%% Find TAD from the result
num = size(Clusters,2);
TAD_Q=[]; % [Number of TD, Average size of TD
count = 1;
for k = Start:End      
    for Ind = 1:2
        Agn = Order_TADNum(Clusters(:,count));
        Assign_Cluster =  Agn';
        %==========================================================================
        %Run the TAD Assessment      
        Extract_TAD_v1; % Extract TAD from clusters identified   
        TAD_Q = [TAD_Q; Q];
        count = count + 1;
    end
    
end
%--------------------------------------------------------------------------
% Write results to file
%--------------------------------------------------------------------------
foldname = [Resultpath,'/Quality'];
if ~exist(foldname, 'dir')
    % Folder does not exist so create it.
    mkdir(foldname);
end
out_path_1 = [foldname,'/'];


TAD_name = strcat(out_path_1,'New_TADQ_',name,'.csv');

dlmwrite(TAD_name,TAD_Q);
