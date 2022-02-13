
Mean_Intra = [];
Mean_Inter = [];
Mean_Pcc = [];
Mean_diff = [];
D = 0;
for k = Start:End   
        for Ind = 1:2
        tdname = [out_path,'KM_',num2str(k),'.',num2str(Ind),'_TAD_BinID.txt'];

    TDname = strcat(tdname); % path to recognized domain
    TD = dlmread(TDname);
    % The average intra-TD and inter-TD contact frequencies
    len = length(TD(:,1)); % find the number of TD
    INTER  = [];
    INTRA = [];
    PCC = [];
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
        % IntraPCC 
        pcc = intraPCC(domain1,Chr_Data);  
        PCC = [PCC;pcc]; % Concatenates the intraPCC results
        %---------------------------------------------

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
    

    %outputa
    fprintf('The average PCC over all TD= %f\n', mean(PCC));
    fprintf('The average IF for over each TD= %f\n', mean(INTRA));
    fprintf('The average_Inter IF for over all TDs= %f\n', mean(INTER));
    fprintf('The average difference= %f\n', mean(out));

    Mean_Intra = [ Mean_Intra ; mean(INTRA)];
    Mean_Inter = [Mean_Inter;mean(INTER)];
    Mean_Pcc = [Mean_Pcc;mean(PCC)];
    Mean_diff = [ Mean_diff;mean(out)];
     if (D < mean(out))
            D =  mean(out);
            Str = ['KM_',num2str(k),'.',num2str(Ind),'_domain.txt'];
     end
        end
end

    intername = [out_path_1 ,'Overall_Inter_Recluster_',Resolution,'.csv'];
    intraname = [out_path_1 ,'Overall_Intra_Recluster_',Resolution,'.csv'];
    diffname = [out_path_1 ,'Overall_Diff_Recluster_',Resolution,'.csv'];
    pccname = [out_path_1 ,'Overall_pcc_Recluster_',Resolution,'.csv'];
    
    dlmwrite(intername,Mean_Inter);
    dlmwrite(intraname, Mean_Intra);  
    dlmwrite(diffname, Mean_diff);
    dlmwrite(pccname, Mean_Pcc);
     %----------------------------------------------------------------------
    % Find the maximum difference
    %----------------------------------------------------------------------
    f= max(Mean_diff);
    readme=strcat(out_path_1,'Readme_Recluster',name,'.txt');
    msg = ['Recommended TAD = ', Str,' with value ',num2str(D)];
    fid = fopen(readme,'wt');
    fprintf(fid, msg');
    fclose(fid);

	if (D >MaxQuality)
       MaxQuality=  D;       
       Outputmsg = ['The name of Best TAD identified Can be found Here:' readme '\n' 'Recommended TAD = ', Str,' with value ',num2str(D) '\n'];
    end



% and the difference between the average intra-TD and inter-TD contact frequencies