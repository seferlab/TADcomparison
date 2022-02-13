% This progrms can find the inter Inter(i, j) denote the
% average of contact frequencies between a bin in TD i and a
% bin in adjacent TD j, where |i–j | = 1
function [Sum,count] = inter_v2(domain_i, domain_j,Chr_Data)

%domain i
Start_i = domain_i(:,1);
End_i = domain_i(:,2);
%domain j
Start_j = domain_j(:,1);
End_j = domain_j(:,2);

Sum = 0;
count = 0 ;
incr = 0;
for i = Start_i:Start_j-1
%     fprintf('Value of i= %d\n',i);
    incr = incr + 1;
    c = 0;
    for j = End_i+1 :  End_j
%     fprintf('Value of i= %d\n',j);
    c = c + 1;
    count =count  + 1;
%     fprintf('%d and %d\n',i,j);
    Sum = Sum + Chr_Data(i,j); 
    if(c==incr)
        break;
    end     
      
    end
end


end



% %domain i
% Start_i = domain_i(:,1);
% End_i = domain_i(:,2);
% %domain j
% Start_j = domain_j(:,1);
% End_j = domain_j(:,2);
% 
% Sum = 0;
% count = 0 ;
% for i = Start_i:Start_j-1
%     for j = End_i+1 :  End_j
%         count =count  + 1;
%         Sum = Sum + Chr_Data(i,j); 
%     end
% end