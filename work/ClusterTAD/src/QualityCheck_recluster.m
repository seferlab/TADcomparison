



    % The average intra-TD and inter-TD contact frequencies
    len = length(TD(:,1)); % find the number of TD
    INTER  = [];
    INTRA = [];
    
    if (len > 1)
        for i = 1:len
            %=============================================
            % For each Domain find intra-TD and inter-TD
            % Intra(i) denote the average of contact frequencies between bins within the same TD i
            % Inter(i, j) denote the% average of contact frequencies between a bin in TD i and a bin in adjacent TD j, where |i–j | = 1
            %---------------------------------------------
            % No 1: Intra    
            domain1 = TD(i,:);     
            Intra = intra(domain1,Chr_Data);
            %=============================================


            % No 2: Inter

            if (i==1)
                 % Adjacent td is just the next one
                domain2 = TD(i+1,:);
               [Sum,count]= inter_v2(domain1,domain2,Chr_Data);
               Inter = Sum/count;
            elseif(i>1 && i < len)
                 % Adjacent td includes prev and next
                domain2 = TD(i-1,:);
                [Sum1,count1]= inter_v2(domain2,domain1,Chr_Data);
                domain2 = TD(i+1,:);
                [Sum2,count2]= inter_v2(domain1,domain2,Chr_Data);
                Inter = (Sum1 + Sum2) / (count1 + count2);
            elseif(i==len)
                 % Adjacent td id just the  prev one     
                domain2 = TD(i-1,:);
                [Sum,count] = inter_v2(domain2,domain1,Chr_Data);
                Inter = Sum/count;
            end
            INTER = [INTER;Inter];
            INTRA = [INTRA;Intra];


        end

        out = INTRA - INTER; % Concaenate the Intra anf Inter Results
    else
        out =  0;
    end
    Intra_Inter_diff  = mean(out);

 