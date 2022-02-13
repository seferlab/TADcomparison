function [ finalD ] = Find_Diagonal(index,T,Limit,zeroRows )
%% This Function Picks the point on the diagonal that falls on a cluster

Index_And_Cluster = [index T]; 
finalD= [];
for i = index(1,1):(index(1,1)+Limit-1)
     if(sum(zeroRows==i)==0)
         D = find(Index_And_Cluster(:,1) == i & Index_And_Cluster(:,2) == i);
         D2 = Index_And_Cluster(D,3);
     else
         D2 = 0;
     end
  
     finalD = [finalD;D2];
end

end

