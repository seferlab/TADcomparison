function mu_s = Mean_density(H,ls,gamma)
% Sub-function of TADs extraction (program: TAD_DP1) via dynamic programming
% Mean scaled density calculation
% H: Input HiC matrix
% ls: size of the blocks
% gamma: Resolution factor: by defaut: 1
% Note that the matrix should be SYMETRIC

% Implemented by
% Jie Chen
% University of Michigan, Ann Arbor
% http://www.jie-chen.com
% dr.jie.chen@ieee.org

%% Default parameter for gamma
if nargin < 3
    gamma =1;
end

L = size(H,1);        % matrix size

%% Mean Density Calculation
% (faster version)
for i = 1 : length(ls)
    ls_gamma = ls(i)^gamma;
    idx = [1:ls(i)];
    Sub_H = H(1:ls(i),1:ls(i));      
    Sub_H = triu(Sub_H)-diag(diag(Sub_H));
    s = sum(Sub_H(:));
    s_all = s / ls_gamma;
    s_pre = s;
    for j = 2 : L-ls(i)+1
        idx_pre = [2:ls(i)]+j-2;
        s =  s_pre - sum(H(j-1,idx_pre)) +  sum(H(j+ls(i)-1,idx_pre));
        s_pre = s;
        s = s / ls_gamma;
        s_all = s_all + s;

    end
    mu_s(i) = s_all/(L-ls(i)+1);
end


% %% Mean Density Calculation
% % (Code according to the expression,slow)
% for i = 1 : length(ls)
%     s_all = 0;
%     for j = 1 : L-ls(i)+1
%         idx = ([1:ls(i)]-1)+j;
%         Sub_H = H(idx,idx);      
%         Sub_H = triu(Sub_H)-diag(diag(Sub_H));
%         s = sum(Sub_H(:));
%         s = s / ls(i)^gamma;
%         s_all = s_all + s;
%     end
%     mu_s(i) = s_all/(L-ls(i)+1);
% end

end