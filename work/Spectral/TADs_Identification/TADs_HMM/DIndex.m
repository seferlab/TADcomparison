function [DI,A,B] = DIndex(H,ls)
% Calculate the directionality index from a Hi-C matrix
% Sub-function for 'TAD_HMM'
% 
%
% Implemented by
% Jie Chen
% University of Michigan, Ann Arbor
% http://www.jie-chen.com
% dr.jie.chen@ieee.org

if nargin < 2
    ls = round(size(H,1)/20);
end


L = size(H,1);


A = zeros(L,1);
B = zeros(L,1);
% downstream
Sub_H = H(1:ls,1:ls);      
Sub_H = triu(Sub_H)-diag(diag(Sub_H));
s = sum(Sub_H(:));
B(1) = s;
for i = 2 : L-ls+1
        idx_pre = [2:ls]+i-2;
        B(i) =  B(i-1) - sum(H(i-1,idx_pre)) +  sum(H(i+ls-1,idx_pre));
end
for i = L-ls+2 : L
    Sub_H = H(i:end,i:end);      
    Sub_H = triu(Sub_H)-diag(diag(Sub_H));
    s = sum(Sub_H(:));
    B(i) = s;
end
B(end)=0;

% upstream
A(ls:end) = B(1:end-ls+1);
for i = 1 : ls-1
    Sub_H = H(1:i,1:i);      
    Sub_H = triu(Sub_H)-diag(diag(Sub_H));
    s = sum(Sub_H(:));
    A(i) = s;
end
A(1) = 0;    
E = (A+B)/2;
DI = (B-A)./abs(B-A).*(((A-E).^2)./E+((B-E).^2)./E);

end