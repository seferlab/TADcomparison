function TAD_mark = TAD_HMM(H,ls)
% TADs extraction via "directionality index" and hidden Markov model
% H: Input HiC matrix
%
% Reference:
% Topological domains in mammalian genomes identified by analysis of chromatin interactions
% J. R. Dixon, S. Selvaraj, et al.
%
% Note: - Personally, I do not think this is a good algorithm. Its modeling
%         was very farfetched.
%       - I recommand the algorithm developped by myself with grapha
%         Laplcian segmentation.
%       - The algorithm developped by D. Filippova via dynamic programming
%         is also not bad, except it would be slow for large scale matrices.
%
% Implemented by
% Jie Chen
% University of Michigan, Ann Arbor
% http://www.jie-chen.com
% dr.jie.chen@ieee.org



%% Directionality index calculation
DI=DIndex(H,ls)';

%% Remove extreme values at two edges
DI(DI>max(DI(ls:end-ls))) = max(DI(ls:end-ls));
DI(DI<min(DI(ls:end-ls))) = min(DI(ls:end-ls));

%% --- Hidden Markov model ---
O = DI;       % set DI as the observations
Q = 3;        % number of states: 3 states.

% Initiliaze initial state
prior0 = [1,0,0]';       % Suppose always starting a TAD

% Initialize the transition matrix
% Initalized intuitively with the transition likelihood
%              1    2     3
transmat0 = [0.1   0.6   0.3       %  1
    0.1   0.6   0.3       %  2
    0.8   0.0   0.2];     %  3

% Loop for number of Gaussian mixtures
for M = 1 : 20
    M
    %     -- Jie's Heuristic initializtion --
    %     mu0 = zeros(size(O,1),M,Q);
    %     mu0(1,:,1) = linspace(max(O), max(O)/20, M);
    %     mu0(1,:,2) =  linspace(max(O)/20, min(O)/20, M);
    %     mu0(1,:,3) =  linspace(min(O)/20, min(O), M);
    %     mu0 = permute(mu0,[1,3,2]);
    %
    %     sp(1) = (max(O)-max(O)/20)/(M-1);
    %     sp(2) = (max(O)/20-min(O)/20)/(M-1);
    %     sp(3) = (min(O)/20-min(O))/(M-1);
    
    %    -- Jie's kmeans initialization --
    mu0 = zeros(size(O,1),M,Q);
    
    sp = zeros(Q,M);
    
    Ot1 = unique(O(O>0));
   if M <= length(Ot1)
        [cls,ct] = kmeans(Ot1, M,'Start','sample');
        mu0(1,:,1) = ct';
        sp(1,:) = kmean_var(Ot1, cls,M);
    else
        mu0(1,:,1) =  linspace(0, min(O), M);
        sp(1,:) = (max(Ot1)/(M-1)/3)^2;
    end
    
    Ot3 = unique(O(O<0));
    if M <= length(Ot3)
        [cls,ct] = kmeans(Ot3, M,'Start','sample');
        mu0(1,:,3) = ct';
        sp(3,:) = kmean_var(Ot3, cls,M);
    else
        mu0(1,:,3) =  linspace(0, min(O), M);
        sp(3,:) = (min(Ot3)/(M-1)/3)^2;
    end
    
    Ot2 = O(O<max(Ot1)/20 & O > min(Ot3)/20);
     if M <= length(Ot2)
        [cls,ct] = kmeans(Ot2, M,'Start','sample');
        mu0(1,:,2) = ct';
        sp(2,:) = kmean_var(Ot2, cls,M);
    else
        mu0(1,:,2) =  linspace(0, min(O), M);
        sp(2,:) = ((max(Ot1)/20 - min(Ot3)/20)/(M-1)/3)^2;
    end
    
    mu0 = permute(mu0,[1,3,2]);
    
    Sigma0 = zeros(size(O,1),size(O,1),M,Q);
    for q = 1 : Q
        for m = 1 : M
            %   Sigma0(1,1,m,q) = sp(q)/2;    % for heuristic initialization
            Sigma0(1,1,m,q) = sp(q,m);    % for kmeans intiaizliaztion
        end
    end
    Sigma0 = permute(Sigma0,[1,2,4,3]);
    
    [ll_trace, prior, transmat, mu, sigma, mixmat,gamma(:,:,M)] = mhmm_em(O, ...
        prior0, transmat0, mu0, Sigma0,1/M*ones(3,M),'max_iter',15);
    
    LL = ll_trace(end);
    AIC(M) = 2*(M*Q + M*Q + Q*(Q-1) + (Q*(M-1)) + (Q-1)) -2*LL;
    % mean  var     A         mixprop     s0
end
[v,idx] = min(AIC);
idx
[v,idx1]=max(gamma(:,:,idx));
figure,
subplot(2,1,1),imagesc(H), colormap hot, colormap(flipud(colormap));
%tbar = zeros(10,size(H,1));
tbar = ones(10,1)*idx1;
subplot(2,1,2),imagesc(tbar)
TAD_mark = idx1;
end


function spv = kmean_var(Op, cls,M)
for i = 1 : M
    spv(i)  = var(Op(cls==i));
end
end