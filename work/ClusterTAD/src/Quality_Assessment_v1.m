%% Find TAD from the result

Len = End - Start + 1;
% Assign_Cluster = zeros(N,Len);
TAD_Q=[]; % Number of TD, Average size of TD
K = Start;
for k = 1:Len
    fprintf ('Result when k = %d\n', K);
    Agn = Order_TADNum(Clusters(:,k));
    Assign_Cluster =  Agn';
    tad_assign=strcat(Resultpath,algorithm,'_',name,'_',num2str(k),'.mat');

    %==========================================================================
    % Quality Assessment
    %--------------------------------------------------------------------------    
    %Run the TAD Assessment      
    Extract_TAD; % Extract TAD from clusters identified   
    TAD_Q = [TAD_Q; Q];
    K = Start + k;
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


TAD_name = strcat(out_path_1,'TADQ_',name,'.csv');
dlmwrite(TAD_name,TAD_Q);
