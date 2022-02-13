function TAD_interval = TAD_DP1(H,gamma)
% TADs extraction via dynamic programming
% H: Input HiC matrix
% gamma: Resolution factor: by defaut: 1

% Reference:
% Multiscale Identification of Topological Domain in Chromatin
% D. Filippova, R. Patro, G. Duggal, C. Kingsford

% Implemented by
% Jie Chen
% University of Michigan, Ann Arbor
% http://www.jie-chen.com
% dr.jie.chen@ieee.org


%% Default parameter for gamma
if nargin < 2
    gamma =1;
end

L = size(H,1);

%% Pre-calculation of mean scaled density
for ls = 2 : L;
   mus(ls) = Mean_density(H,ls,gamma);
end

%%
y_opt(1) = 0;
for l = 3 : L
    for k = 2: l-1
        y_opt_cand(k-1) = y_opt(k-1) + max(Scaled_density(H,k,l,gamma) - mus(l-k),0);
    end
    [y_opt(l),pos0(l)] = max(y_opt_cand);
end

c = 1;
pos_rc(c) = pos0(L);
while 1
   if pos_rc(c)<=2  |  pos_rc(c) == pos0(pos_rc(c)) 
       break;
   else
       c = c+1;
       pos_rc(c) = pos0(pos_rc(c-1));
   end

end 
TAD_interval = fliplr(pos_rc);
end