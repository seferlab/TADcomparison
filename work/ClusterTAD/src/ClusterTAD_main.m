%%========================================================================
 % Implementation of ClusterTAD algorithm
 %  Language: MATLAB
	
 %	@author: Oluwatosin Oluwadare 
 %	Department of Electrical Engineering and Computer Science
 %	University of Missouri, Columbia
 %	Columbia
 %	USA	
 %	email: 	oeow39@mail.missouri.edu
	
 %	Last Update: 2/7/2018
	
 %	Instructions: https://github.com/BDM-Lab/ClusterTAD



%clear;
warning('off','all');
%  =============================================================================================
%   ClusterTAD  Variables
%  =============================================================================================

%	change variable "filepath" to the Input data file path	
%	change variable "name" to the Input data name
%	change variable "chromo" to the chromosome data name
%	change variable "Reso" to the Input data Resolution		
%	change variable "Option" to specify  the  K finding method. if option == 0, use rough estimate and if option == 1, use elbow method (Read about this in our manuscript)
%	specify the maximum TAD size : change variable "Max_TADsize" (optional)

%  ============================================================================================
MaxQuality = 0;			   % Get the Maximum TAD Quality

%filepath='data\';          % filepath

%name= '30x30_data.txt';          % filename

%chromo = 'chr';             % chromosome name

%Res = 40000;                  % 40000 = 40KB , 100000 = 100KB

KB = 1000;                    % KB = 1000 :: Constant

Option = 0;                % Option: if option == 0, use rough estimate and if option == 1, use elbow method

%Max_TADsize = 800000;          % Maximum TAD size: 800KB

%outputfolder_name = 'data_Results';

algorithm = 'KM';          % algorithm used code: KM == Kmeans Algorithm
%  ============================================================================================
Resolution=[num2str(Res/KB),'kb'];
filename = strcat(filepath,name);
Chr_Data = dlmread(filename); 

fprintf('Data set loaded.\n');

[path,name,ext] = fileparts(name);
%==========================================================================
% Make directory if it doesn't exist
%--------------------------------------------------------------------------
foldname = [outputfolder_name,'/',Resolution,'/',chromo];
if ~exist(foldname, 'dir')
    % Folder does not exist so create it.
    mkdir(foldname);
end

Resultpath = [foldname,'/'];
N = length(Chr_Data);
% Rough estimate
M = floor(sqrt(N/2));
K = [];
% define the number of cliuster[ Using Elbow method ]
for i = 1: (M + 20)
    K = [K,i];
end

% Get dimension
d = length(K);
distortion=zeros(d,1);
test_num = 4; %Perform clustering about 5 times ans select minimum

%  =========================================
% Find Zero Rows
%  =========================================
zeroRows = [];
for i = 1:size(Chr_Data,1)
    if isequal(Chr_Data(i,:),zeros(1,size(Chr_Data,2)))
        zeroRows = [zeroRows; i];    
    end
end
fprintf('The number of Zero Rows = %d \n', length(zeroRows));

% ============================================
% Create format
% ============================================
Limit = length(Chr_Data);
inLimit = 1;
index = [];
Comb_Data=[];
for i  = 1:Limit    
    count = 0;  
    for j = i  
        if(sum(zeroRows==i)==0)
             count = count +1 ;% keep a count to use only the useful dataset
             if(count<=inLimit)
                 D  = [Chr_Data(i,1:length(Chr_Data)) Chr_Data(1:length(Chr_Data),j)' ];         
                 Comb_Data = [Comb_Data ; D];
                 in= [i j];
                 index = [index ;in];
             end
        end
       
    end         
end  

Comb_Data(~isfinite(Comb_Data))=0;
% ============================================
% Decide on the K value to use for clustering 
% ============================================


if ( Option ==0)
    k_opt =  M;
    Start = k_opt- 1;
    End = k_opt + 1;
    
    if (Start < 1 )
        Start = 2;
    end
     if (Start > Limit )
        End = k_opt + 2;
    end
     row = size(Comb_Data,1); leng = End - Start + 1;
    Clusters = zeros(row,leng);
    fprintf('kmeans_Clustering started........\n'); 
    K = Start;
    for i = 1:leng
        No_of_Cluster = K;   % No of clusters        
        fprintf('kmeans_Clustering :: K = %d ........\n', No_of_Cluster);    
        [T,~,sumd] = kmeans(Comb_Data,No_of_Cluster,'Distance','sqEuclidean','start','uniform', 'emptyaction',...
            'drop'); 
        % try differnet tests to find minimun disortion under k_temp clusters
        destortion_temp=sum(sumd);
        Clusters(:,i) = T;
        for test_count=2:test_num     
            [V,~,sumd]=kmeans(Comb_Data, No_of_Cluster,'Distance','sqEuclidean','emptyaction','drop');
            if ( sum(sumd) < destortion_temp )
               Clusters(:,i) = V;  
               destortion_temp = sum(sumd);
            end
        end
        K = K + 1;
        

    end
    
else

    % Clustering using Elbow : implementing kmeans
    row = size(Comb_Data,1);
    C = zeros(row,length(K));
    fprintf('kmeans_Clustering started........\n');    
    for i = 1:length(K)
        No_of_Cluster = K(i);   % No of clusters    
        fprintf('kmeans_Clustering :: K = %d ........\n',i);    
        [T,~,sumd] = kmeans(Comb_Data,No_of_Cluster,'Distance','sqEuclidean','start','uniform', 'emptyaction',...
            'drop'); 
        % try differnet tests to find minimun disortion under k_temp clusters
        destortion_temp=sum(sumd);
        C(:,i) = T;
        for test_count=2:test_num     
            [V,~,sumd]=kmeans(Comb_Data, No_of_Cluster,'Distance','sqEuclidean','emptyaction','drop');
            if ( sum(sumd) < destortion_temp )
               C(:,i) = V;  
               destortion_temp = sum(sumd);
            end
        end
        distortion(i,1)=destortion_temp;  

    end


    variance=distortion(1:end-1)-distortion(2:end);
    distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
    [r,~]=find(distortion_percent>0.8);
    k_opt=r(1,1)+1;

    Start = k_opt- 5;
    End = k_opt + 5;
     if (Start < 1 )
        Start = 2;
    end
     if (Start > Limit )
        End = k_opt + 2;
    end
    row = size(Comb_Data,1); leng = End - Start + 1;
    Clusters = zeros(row,leng);
    ind = Start;
    for i = 1:leng        
      Clusters(:,i) = C(:,ind);
      ind = ind + 1;
    end
    
    dlmwrite([Resultpath,'kmeans_distortion_percent.txt'],distortion_percent);
    dlmwrite([Resultpath,'kmeans_variance.txt'],variance);
    dlmwrite([Resultpath,'kmeans_distortion.txt'],distortion);


    %% Obtaining the second derivative of the distortion
    derivative = diff(distortion,2);

    figure ; 
    plot(distortion_percent,'b*--');
    xlabel('Number of Clusters(K)') % x-axis label
    ylabel('Percentage of variance (\%)') % y-axis label
    title('Elbow Method Plot1');

    figure ;
    plot(distortion);
    xlabel('Number of Clusters(K)') % x-axis label
    ylabel('Average within Cluster Distance to Centoid') % y-axis label
    title('Elbow Method Plot2');
    

end
fprintf('Recommended Number of Cluster (K) = %d \n',k_opt);

readme=strcat(Resultpath,'Readme_',name,'.txt');
msg = ['Recommended Number of Cluster (K) = ', num2str(k_opt)];
fid = fopen(readme,'wt');
fprintf(fid, msg');
fclose(fid);

%-------------------------------------------------------------------------
% Save Cluster assignment Alone
foldname = [Resultpath,'/Clusters'];
if ~exist(foldname, 'dir')
    % Folder does not exist so create it.
    mkdir(foldname);
end
out_path_1 = [foldname,'/'];
clustfile=strcat(out_path_1,algorithm,'_Cluster_',name,'_.txt');

%-------------------------------------------------------------------------
% Label Clusters

Cluster = zeros(Limit,leng);
for j = 1:leng
     for i = 1:length(index)
          Cluster(index(i,1),j) = Clusters(i,j);
     end
end


dlmwrite(clustfile,Cluster);
Clusters = Cluster;
%--------------------------------------------------------------------------
disp('======================== Quality Assessment==========================');

%% Perform Quality Assessment
Quality_Assessment_v1;
QualityCheck;
disp('===================== Quality Assessment Completed==========================');

disp('===================== Re-Cluster Points=====================================');
Recluster_v1;

%  Optional
% Note: Experiments showed that the Quality of TADs produced at 
%       2nd Reclustering (Recluster_v2) equivalent to 3rd Clustering operation sometimes doesn't improve performance
% Recluster_v2;  
disp('===================== Re-Clustering Completed =============================');
disp('find the results of the Quality Assessment in the Quality/ directory');
fprintf( Outputmsg);
disp('Find the TADs identified in the TAD/ directory');
disp('===================== ClusterTAD Completed =============================');
final_file=[outputfolder_name,'/',Resolution,'/',chromo,'/TADs/',Str];
save_file=[outputfolder_name,'/',Resolution,'/',chromo,'/final_domain.txt'];
movefile(final_file,save_file);