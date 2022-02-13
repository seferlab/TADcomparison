function s = Scaled_density(H,k,l,gamma)
% Sub-function of TADs extraction (program: TAD_DP1) via dynamic programming
% Scaled density calculation
% H: Input HiC matrix
% k: starting loci
% l: ending loci
% gamma: Resolution factor: by defaut: 1

% Implemented by
% Jie Chen
% University of Michigan, Ann Arbor
% http://www.jie-chen.com
% dr.jie.chen@ieee.org

%% Default parameter for gamma
if nargin < 4
    gamma =1;
end

%% Density Calculation
Sub_H = H(k:l,k:l);
Sub_H = triu(Sub_H)-diag(diag(Sub_H));
s = sum(Sub_H(:));
s = s / (l-k)^gamma;

end