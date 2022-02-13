diary off;
nn=  ['_',num2str(k),'__',num2str(Ind)];
foldname_1 = [Resultpath,'/TADs'];
if ~exist(foldname_1, 'dir')
    % Folder does not exist so create it.
    mkdir(foldname_1);
end
out_path = [foldname_1,'/'];
filename = [out_path,algorithm,nn,'TADReport_.txt'];
diary (filename)
diary on;
fprintf('TAD for %s CLUSTERING\n',algorithm);
figure;
n = size(Chr_Data,1); % specify the size of matrix to display :::default
% n = 100;
Visualize(log(Chr_Data(1:n,1:n)));
newB = F_TAD(Chr_Data(1:n,1:n), Assign_Cluster, algorithm, nn,out_path ,Res );
avg_size = (sum(newB(:,2))- sum(newB(:,1)))/size(newB,1);
fprintf('Average size = %f\n',avg_size);
diary off;
%close;
Q = [length(newB(:,1)) avg_size]; % Number of TD, Average size of TD
title_text = sprintf('TD for %s Implementation',algorithm);
title(title_text)
