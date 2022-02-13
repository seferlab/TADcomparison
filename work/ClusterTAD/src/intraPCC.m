function [ Intra_Average] = intraPCC( domain,Chr_Data )
%% Pearson’s correlation coefficient (PCC) between the contact profiles of bins within a TD as a quality measurement
Start = domain(:,1);
End = domain(:,2);
Mat = [];
for i = Start:End
    temp = [];
    for j = Start : End
        temp = [temp Chr_Data(i,j)];       
    end
    Mat = [Mat;temp];
    
end
% Now that Matrix has been formed, do a correlation
Matrix = corr(Mat); % Matrix contains the pearsons correlation of each bin element
Matrix(isnan(Matrix)) = -0.2;
%-----------------------------------------------------------------------------
% The average intra-TD Pearson’s correlation coefficient (PCC) and weighted PCC
%--------------------------------------------------------------------------
%1. Intra PCC
 num = size(Mat,2);
 Intra_Average = sum(sum(triu(Matrix)))/((num+1)*num/2); % Average

end

