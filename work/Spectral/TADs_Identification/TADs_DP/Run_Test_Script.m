% Script for testing topological domain extraction algorithm. 
% Call the dynamic programming method.
%
% Implemented by
% Jie Chen
% http://www.jie-chen.com
% dr.jie.chen@ieee.org

clear, close all, clc

% Algorithm selection Alg = 1 or 2
Alg = 1;
% Fiedler value threshold
gamma = 0.35;

%% Load the example data
load ../HiC_chr22_100kb

%% Remove unmappable region
idx = sum(H)~=0;
H = H(idx,idx);

%% Apply a transformation
Ht = ceil(H);
% log transformation, and saturated by 6;
Ht = min(log(Ht),6);
% Process -inf, because log(0) = -inf
Ht(Ht == -inf) = -1;
% Shift to be positive
Ht = Ht + 1.001;

% or we can use:
%Ht = sqrt(sqrt(H));


%% Call Algorithm
TAD_boundaries = TAD_DP1(Ht,gamma);


% Display
Draw_TADs(Ht, TAD_boundaries,[0,6]);
