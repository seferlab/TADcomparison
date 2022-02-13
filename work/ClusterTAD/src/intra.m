% This find the intra-TDcontact frequency and returns average
function [Intra_Average,Sum,count] = intra(domain,Chr_Data)
Start = domain(:,1);
End = domain(:,2);
Sum = 0;
count = 0 ;
for i = Start:End
    for j = i+1 : End
       count =count  + 1;
       Sum = Sum + Chr_Data(i,j); 
    end
    
end

% Find the average

if (Sum >  0 && count > 0)
    Intra_Average = Sum/count;
else
    Intra_Average = 0;
end

end