function sumVec=sumDiag(A,antiD)
%sumVec=sumDiag(A,antiD)
%Summation across all separate diagonals (or antidiagonals) of matrix A returned in 
%a column vector. Works for 2D matrices of any size. Does not involve any
%for-loops. The first element of the returned vector is the upper right element 
%(if diagonals) or upper left element (if antidiagonals), and so on.
%
%Input:
% A       -  2D Matrix
% antiD   -  Optional. Set = 1 to sum antidiagonals. Default = 0 (diagonals).
%Output:
% sumVec  -  Vector of diagonal sums
%
%Written by: Marcus Björk, Uppsala University, 2013

[N,M]=size(A);
endFlip=0;

%Alternative to sum antidiagonals (default - diagonals)
if nargin==1
    antiD=0;
end

%Make sure matrix is not wide (reduces computaional burden for very wide matrices)
if M>N
    A=A.';
    [N,M]=size(A);
    endFlip=1;
end

Amod=zeros(N+M-1,M);
logVec = [false(M-1,1);true(N,1);false(M-1,1)];
indMat = bsxfun(@plus, (1:M+N-1)',0:M-1);
if antiD
    indMat=flipud(indMat);
    %indMat = bsxfun(@plus, (M+N-1:-1:1)',0:M-1);
end
logMat = logVec(indMat);    
Amod(logMat)=A;

sumVec=sum(Amod,2); %Perform summation

%Return in correct order (compensate for wide matrix being transposed)
if endFlip
    sumVec=flipud(sumVec);
end
