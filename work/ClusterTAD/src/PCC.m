function [ output_args ] = PCC()
%% Find the PCC quality check 2
name = 'nij.chr17';
tdname = 'chr17_new.txt';
%% Quality Comparism
% This function does the quality comparism of each TAD found by the Different
% methods using Approach Described in TopDom
% Intra-Inter Difference Comparisom
% Intra PCC Quality
% WPCC

%% Load the normalized Contact Matrix
filepath = 'C:\Users\Tosin\Documents\Research\TAD\Dataset\hESC.norm\hES\nij\';
% filepath = 'C:\Users\Tosin\Dropbox\Research\TAD\Quality\';
path ='C:\Users\Tosin\Dropbox\Research\TAD\Quality\';
filename = strcat(filepath,name); %filepath is path to chromosome data
Chr_Data = dlmread(filename);

% Load the TD
TDname = strcat(path, tdname); % path to recognized domain
TD = dlmread(TDname);

% The average intra-TD and inter-TD contact frequencies
len = length(TD(:,1)); % find the number of TD

INTRA = [];
for i = 1:len
    %=============================================
    % IntraPCC 
    domain = TD(i,:);   

    Intra = intraPCC(domain,Chr_Data);  
    INTRA = [INTRA;Intra]; % Concatenates the intraPCC results
    %=============================================

end

output = 'C:\Users\Tosin\Dropbox\Research\TAD\Quality\CHDF\';
Out_name = strcat(output,name,'_intraPCC');

%dlmwrite(Out_name, INTRA);

fprintf('The average = %d\n', mean(INTRA));

end





% the intra-TD Pearson’s correlation coefficient (PCC) 
% the intra-TD weighted Pearson’s correlation coefficient (wPCC),
% and the difference between the average intra-TD and inter-TD contact frequencies