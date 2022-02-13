function [newB] = F_TAD_recluster(Chr_Data, Assign_Cluster,base_ind)
count = 1;
Border = [];
start = 1;
Limit = length(Assign_Cluster);
for i = 2:Limit
    if(Assign_Cluster(i) ~= Assign_Cluster(start))
        Border = [Border;[start,i-1]];
        start = i;
        count = count + 1;
    end
end
Border = [Border;[start,Limit]];

 newB = [];
    for j = 1:length(Border(:,1))
        if ((Border(j,2)- Border(j,1) + 1) > 3)
          newB = [ newB ; Border(j,:) ] ;    
        end
    end

newB = newB  + base_ind;
%===========================================================================
% Find the gaps classified as domain
zeroRows = [];
for i = 1:size(Chr_Data,1)
    if isequal(Chr_Data(i,:),zeros(1,size(Chr_Data,2)))
        zeroRows = [zeroRows; i];    
    end
end
%===========================================================================
Redefine = [];
newB_domain = [];
B_end = newB(:,2);
for i = 1:length(B_end)
 if(sum(zeroRows==B_end(i))==0)
     X = [newB(i,:)  1];   % Recognized as domain = 1
     newB_domain = [newB_domain ; newB(i,:)];
 else
     X = [newB(i,:)  0];   %  Recognized as void  = 0
 end
     Redefine = [Redefine; X];
end

newB = newB_domain;

fprintf('The Number of TAD in re-cluster= %d\n', length(newB(:,1)));

end

