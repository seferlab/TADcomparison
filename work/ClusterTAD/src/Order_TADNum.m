function [ Assign]  = Order_TADNum(Found_TAD )
%% This Function number the TAD clusters found.
%--------------------------------------------------------------------------
len = length(Found_TAD); % Length of the current TAD
count = 1;
Assign(1) = count; 
for i = 2:len
    if(Found_TAD(i-1) == Found_TAD(i))
        count = count + 0;
    else
        count = count + 1;
    end
        Assign(i) = count;
end

end

