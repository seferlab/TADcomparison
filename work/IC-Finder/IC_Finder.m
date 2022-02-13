function [dom, Z, pd, pb] = IC_Finder(HiCmap, varargin)
%
%   by Noelle HADDAD, Cédric VAILLANT and Daniel JOST
%   May 2016
%   Please read the LICENSE file before using IC_Finder
%   Any remarks and suggestions should be addressed to Daniel JOST daniel.jost@imag.fr
%
%IC_Finder allows to segment HiC maps into interacting compartments (IC).
%   Required toolbox : Statistics and Machine Learning Toolbox for MATLAB users
%                      statistics for OCTAVE users (pkg load statistics)
%
%   DOM = IC_FINDER(HICMAP), where HICMAP is a symmetric square matrix or 
%   the path to a .txt or .mat file containing it. It creates a two-columns
%   matrix DOM defining the boundaries of detected domains. The
%   segmentation is based on two intuitive parameters sigmaMinus and 
%   sigmaPlus whose default values (resp. 0.3 and 3) have been learned to 
%   give optimal results on a large variety of experimental HiC data. To
%   change sigmaMinus and/or sigmaPlus, see below.
%   In DOM, each line represents one domain, in the first and second
%   columns, there are respectively the first and last bins of the domain.
%   DOM is by default saved as a text file in the same directory of HICMAP.
%   The file name is by default the same as HICMAP with the extension :
%   ' _domains.txt ' if HICMAP is a path or is card_date_time_domains.txt
%   if HICMAP is a matrix. To change the path and/or the name of DOM or to
%   not save DOM as a text file see options below.
%   Note that if from a line to another there are missing bins in DOM,
%   it means that for these bins there is not enough signal to segment them.
%   Example for DOM:
%   1 10
%   11 15
%   20 32
%   From bins 16 to 19 include, the segmentation is not possible. 
%
%   IC_FINDER ignores columns and rows of HICMAP that contain more than 75%
%   missing values (0 or NaNs) on the 20 nearest coefficients of the
%   diagonal.
%
%   [DOM, Z] = IC_Finder(HICMAP) returns in addition a matrix Z that encodes
%   for a tree of hierarchical clusters. By using
%   the function 'dendrogram', it is possible to plot the 
%   hierarchical binary cluster tree. Syntax : dendrogram(Z)
%
%   [DOM, Z, PD, PB] = IC_Finder(HICMAP, 'Option', 'resampling')
%   performs statistical resampling of the input HiC map to estimate PD,
%   the probability that two bins belong to the same domain and PB,
%   the probability for a bin to be predicted as a IC boundary.
%   PD is a matrix which has the size of HICMAP. PD(i,j) is the probability
%   that bins i and j colocalize in the same domain. PB is a 2 columns matrix.
%   In each line, the first column corresponds to the number of the bin,
%   the second column to the boundary probability.
%   Matrices PD and PB are saved in text files by default, to not
%   save them set 'SavePD' and 'SavePB' to 0 (see below).
%   The resampling is done Nr times with a noise level defined by the
%   parameter Nl. The higher is Nl the higher is the noise on the resampled
%   maps. By default, Nr=100 and Nl=1. These 2 values can be modified (see 
%   below).
%
%   [DOM, Z] = IC_Finder(HICMAP, 'Option','hierarchy') performs
%   segmentation based on SigmaZero, an optional parameter controlling the 
%   merging of two clusters. By default, SigmaZero=5 if 'Option'='hierarchy'
%   but can be tuned by the user (see below).
%
%   [DOM, Z, PD, PB] = IC_Finder(HICMAP, 'Option','rh') with 'rh' the
%   abbreviation for 'resampling and hierarchy'. It couples the options
%   'resampling' and 'hierarchy' (see above).
%
%   [ ... ] = IC_Finder(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to quantify the reliability of the 
%   predicted domains (resampling) and/or to infer higher-order levels
%   of chromatin organization (hierarchy) and/or to save the data in text
%   files. Parameters and values are not case-sensitive. Parameters are : 
%
%   'Option'  -  Choices are :
%      'no'           - no option (the default)
%      'resampling'   - statistical resampling of HICMAP (see above)
%      'hierarchy'    - inference of higher-order levels of chromatin
%                       organization (see above)
%      'rh'           - Resampling and Hierarchy (see above)
%
%   'Nr'  -  Number of resampled maps. A positive integer between 10 and 
%            500, default is 100 if 'Option'='resampling' or 'Option'='rh'.
%            Otherwise, default Nr is equal to 1.
%
%   'Nl'  -  Noise level to resample HICMAP. A float between 0.1 and 2.
%            Default is 1 to represent the typical number of contacts in 
%            recent HiC experiments.
%
%   'SigmaZero'  -  Positive float to control the merging of two clusters
%            in the 'hierarchy' or 'rh' options. By default, SigmaZero=5.
%            If SigmaZero is a vector, the segmentation will be
%            done for all  number specified in SigmaZero. By this way, one
%            can move through hierarchy of the domains. In this
%            case, outputs DOM, Z, PD and PB have a third dimension which 
%            is equal to length(SigmaZero). The third dimension allows to 
%            save DOM, Z, PD and PB for each sigmaZero value.
%              
%   Note that if you omit to specify 'Option' but you specify Nr, Nl or 
%   SigmaZero a warning message indicating that 'Option' has been  
%   automatically changed will be printed. 
%               
%   'SigmaMinus' and 'SigmaPlus'  -   Low and high thresholds to control 
%            the merging of two clusters in the default mode or in option 'resampling'.
%            SigmaMinus and SigmaPlus are two independent parameters.
%            By default, SigmaMinus=0.3 and SigmaPlus=3.
%                   
%   'SaveDomains' - 1 (default) or 0 to respectively save domains in a file
%            text or not. The filename is 'XX_domains.txt' where XX is the
%            parameter 'Name'. (See below the default value for 'Name')
%                   
%   'SaveZ'       - 0 (default) or 1 to save Z in a file text or not. The 
%            filename is 'XX_Z.txt' where XX is the the parameter 'Name'.
%
%   'SavePD       - 1 (default) or 0 to save PD in a file text or not. The
%            filename is 'XX_PD.txt' where XX is the the parameter 'Name'.
%
%   'SavePB       - 1 (default) or 0 to save PB in a file text or not. The
%            filename is 'XX_PB.txt' where XX is the the parameter 'Name'.
%
%   'PlotFigures  - 1 (default) or 0 to plot results in figures. Fig1 shows
%            the HiCmap with found boundaries. Fig2 adds the dendrogram 
%            representing the hierarchical clustering. Fig3 presents PB and
%            PD if option "Resampling" is valid. These 3 figures are
%            plotted for each level of hierarchy (ie. for each value in 
%            SigmaZero if option "Hierarchy" is valid).
%
%   'SaveFigures  - 1 (default) or 0 to save plotted figures (pdf format).
%
%   'Path'         - Path to save outputs. By default, path is the one of 
%            HICMAP if HICMAP is a string. By default, path is the current 
%            directory if HICMAP is directly a matrix. 
%
%   'Name'         - Name to save outputs. By default, if HICMAP is a 
%            string, Name is equal to the HICMAP filename without the extension.
%            If HICMAP is not a path but a matrix Name is by default :
%            card_YearMonthDay_Hour:Minute:Second .
%
%   Examples:
%
%      % To detect domains on maptest_1.txt and save maptest_1_domains.txt
%      ( in example/ folder )
%      dom = IC_Finder('example/maptest_1.txt');
%      or to not save the result
%      card=load('example/maptest_1.txt');
%      dom = IC_Finder(card,'SaveDomains',0,'SaveFigures',0);
%
%      % To detect domains using an optional threshold and to save results in
%      maptest_1_domains_sigmaZero6.txt
%      dom = IC_Finder('example/maptest_1.txt','Option','hierarchy','SigmaZero',6);
%      or
%      IC_Finder('example/maptest_1.txt','SigmaZero',6);
%      or
%      IC_Finder('example/maptest_1.txt','SigmaZero',6,'path','example/');
%
%      % To detect domains and to compute PD & PB with default parameters. It
%      creates 4 files in example/ folder : maptest_1_domains_nr100_nl1.txt
%      maptest_1_Z_nr100_nl1.txt maptest_1_pd_nr100_nl1.txt 
%      maptest_1_pb_nr100_nl1.txt
%      [dom,Z,pd,pb] = IC_Finder('example/maptest_1.txt','Option','resampling','SaveZ',1);
%
%   See also PIC_FINDER, DENDROGRAM.
%


% Inputs
nameD = 0;
if ischar(HiCmap)
    C = load(HiCmap); 
    % path and name without extension to save outputs at the end
    [pathstr,defaultName] = fileparts(HiCmap);
    if ~isempty(pathstr)
        defaultPath=[pathstr '/'];
    else
        defaultPath=pathstr;
    end
else
    C = HiCmap;
    defaultPath='';
    current_date = clock;
    name=['card_' num2str(current_date(1)) num2str(current_date(2),'%02.0f') num2str(current_date(3),'%02.0f')...
        '_' num2str(current_date(4),'%02.0f') ':' num2str(current_date(5),'%02.0f') ':' num2str(current_date(6),'%02.0f')];
    defaultName=strrep(name,' ','_');
    nameD = 1; % if name is date
end

% if sparse format for C
if size(C,2)==3 && size(C,1)>3
    Cs=C; % Csparse
    N=max(Cs(:,1:2));
    C=zeros(N);
    for i=1:size(Cs,1)
        C(Cs(i,1),Cs(i,2))=Cs(i,3);
        C(Cs(i,2),Cs(i,1))=Cs(i,3);
    end
end

p = inputParser;
p.FunctionName = 'IC_Finder';

defaultOption = 'no';
validOption = {'no','resampling','hierarchy','rh'};
checkOption = @(x) any(validatestring(x,validOption));
defaultNr = 100; 
defaultNl = 1; 
defaultSigmaZero = 5;
defaultSigmaMinus = 0.3; % Default parameters for threshold
defaultSigmaPlus = 3;
defaultSaveDomains=1; 
defaultSaveZ=0; 
defaultSavePD=1; 
defaultSavePB=1;
defaultPlotFigures=1;
defaultSaveFigures=1;
nc = 904;

addRequired(p,'HiCmap');
addParamValue(p,'option',defaultOption,checkOption)
addParamValue(p,'nr',defaultNr,@isnumeric)
addParamValue(p,'nl',defaultNl,@isnumeric) 
addParamValue(p,'sigmaZero',defaultSigmaZero,@isnumeric)
addParamValue(p,'sigmaMinus',defaultSigmaMinus,@isnumeric)
addParamValue(p,'sigmaPlus',defaultSigmaPlus,@isnumeric)
addParamValue(p,'saveDomains',defaultSaveDomains)
addParamValue(p,'saveZ',defaultSaveZ)
addParamValue(p,'savePD',defaultSavePD)
addParamValue(p,'savePB',defaultSavePB)
addParamValue(p,'plotFigures',defaultPlotFigures)
addParamValue(p,'saveFigures',defaultSaveFigures)
addParamValue(p,'path',defaultPath)
addParamValue(p,'name',defaultName)


parse(p,HiCmap,varargin{:})

nr = p.Results.nr;
nl = p.Results.nl; 
sigmaZero = p.Results.sigmaZero;
sigmaMinus = p.Results.sigmaMinus;
sigmaPlus = p.Results.sigmaPlus;
option = p.Results.option;
saveDomains = p.Results.saveDomains;
saveZ = p.Results.saveZ;
savePD = p.Results.savePD;
savePB = p.Results.savePB;
plotFigures = p.Results.plotFigures;
saveFigures = p.Results.saveFigures;
Path = p.Results.path;
name = p.Results.name;

def = p.UsingDefaults;

% Inputs & warnings
validateattributes(nr,{'numeric'},{'positive'},'IC_Finder','Nr')
validateattributes(nl,{'numeric'},{'positive'},'IC_Finder','Nl')
validateattributes(sigmaZero,{'numeric'},{'nonnegative'},'IC_Finder','SigmaZero')
validateattributes(sigmaMinus,{'numeric'},{'nonnegative'},'IC_Finder','SigmaMinus')
validateattributes(sigmaPlus,{'numeric'},{'nonnegative','>=',sigmaMinus},'IC_Finder','SigmaPlus')
validateattributes(saveDomains,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','SaveDomains')
validateattributes(saveZ,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','SaveZ')
validateattributes(savePD,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','SavePD')
validateattributes(savePB,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','SavePB')
validateattributes(plotFigures,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','PlotFigures')
validateattributes(saveFigures,{'numeric'},{'scalar','>=',0,'<=',1},'IC_Finder','SaveFigures')
validateattributes(Path,{'char'},{})
validateattributes(name,{'char'},{})

% Nothing specified
if ( sum(strcmp(def,'option'))==1 && sum(strcmp(def,'nr'))==1 ...
        && sum(strcmp(def,'nl'))==1 && sum(strcmp(def,'sigmaZero'))==1 )
    nr=1; sigmaZero=0;
% Only 'Option' specified (no Nr, Nl and SigmaZero)
elseif ( sum(strcmp(def,'option'))~=1 && sum(strcmp(def,'nr'))==1 ...
        && sum(strcmp(def,'nl'))==1 && sum(strcmp(def,'sigmaZero'))==1 )
    if sum(strfind(option,'no'))
        nr=1;  sigmaZero=0;
    elseif sum(strfind(option,'re'))
        sigmaZero=0;
    elseif sum(strfind(option,'hi'))
        nr=1; 
    elseif sum(strfind(option,'rh'))
    end
% Nr and/or Nl specified    
elseif ( sum(strcmp(def,'nr'))~=1 || sum(strcmp(def,'nl'))~=1 ) && sum(strcmp(def,'sigmaZero'))==1
    if sum(strfind(option,'no'))
        warning('Option is put to ''resampling'' because you specified Nr or Nl.');
        sigmaZero = 0;
    elseif sum(strfind(option,'re'))
        sigmaZero = 0;
    elseif sum(strfind(option,'hie'))
        warning('Option is put to ''rh'' because you specified Nr or Nl.');
    end
    if nr>500 || nr<10  
        warning('Nr should be an integer in [10, 100]');
    end
    if nl<0.1 || nl>2
        warning('Nl should be in [0.1, 2]');
    end
% SigmaZero specified
elseif (sum(strcmp(def,'nr'))==1 && sum(strcmp(def,'nl'))==1 ) && sum(strcmp(def,'sigmaZero'))~=1
    if sum(strfind(option,'no'))
        warning('Option is put to ''hierarchy'' because you specified SigmaZero.');
        nr = 1;
    elseif sum(strfind(option,'res'))
        warning('Option is put to ''rh'' because you specified SigmaZero.');   
    elseif sum(strfind(option,'hi'))
        nr = 1; 
    end
% (Nr or Nl ) and SigmaZero specified 
elseif ( sum(strcmp(def,'nr'))~=1 || sum(strcmp(def,'nl'))~=1 ) && sum(strcmp(def,'sigmaZero'))~=1
    if sum(strfind(option,'no')) || sum(strfind(option,'hi')) || sum(strfind(option,'res'))
        warning('Option is put to ''rh'' because you specified Nr or Nl and SigmaZero.');
    end
    if nr>500 || nr<10  
        warning('Nr should be an integer in [10, 100]');
    end
    if nl<0.1 || nl>2
        warning('Nl should be in [0.1, 2]');
    end
end

% Warning on filename if not specified
if nameD==1 && ( saveDomains==1 || saveZ || (savePD && nr>1) || (savePB && nr>1) ) && sum(strcmp(def,'name'))==1
warning(['As you did not specified ''Name'' nor path for HiCmap in inputs, outputs are going to be saved in the current folder with the name :  '''...
    name '''. If you do not want to save outputs see description of IC_Finder.m']);
end

% Put 0 to the diagonal of the contact map
Cinitial=C; Ninitial = size(Cinitial,1);
for i=1:Ninitial
    C(i,i) = 0;
end

% Put 0 instead of nan
C(isnan(C))=0;

% Delete rows or columns with more than 75% of zeros around diagonal
notEnoughSignal = [];
for i=1:Ninitial
    st=max(1,i-10);
    en=min(Ninitial,i+10);
    if sum(C(st:en,i)>0)/length(st:en)<0.75
        notEnoughSignal(end+1)=i;
    end
end
C(notEnoughSignal,:)=[];
C(:,notEnoughSignal)=[];
        
N = size(C,1);

alpha = 0;
for i=1:N-1
    alpha=alpha+C(i,i+1);
end
alpha=alpha/(N-1);
f=nl*nc/alpha;

Zoutput = zeros(Ninitial-1,3,length(sigmaZero));
domoutput = zeros(1,2,length(sigmaZero));
pdoutput = zeros(Ninitial,Ninitial,length(sigmaZero));
pboutput = zeros(1,2,length(sigmaZero));

for ith=1:length(sigmaZero)
    
    b = [];
    pd = zeros(Ninitial);
    
    for ir=1:nr
        
        if ir==1
            card0 = C;
        else
            card0=f*C;
            card0=poissrnd(card0);
            for ii=1:N
                card0(ii,ii)=0;
                for jj=ii+1:N
                    card0(ii,jj)=card0(jj,ii);
                end
            end
        end
        
        % Normalisation factor card (useful to compute dt) & coeff for weighted
        % mean linkage
        cm=zeros(N-1,1);
        for i=1:N-1
            d=diag(card0,i);a=find(d);d=d(a);a=isfinite(d);a=find(a);
            cm(i)=mean(d(a));
        end
        
        coeff=zeros(N); % used with wm linkage
        cb=zeros(N);
        for i=1:N
            for j=i+1:N
                cb(i,j)=card0(i,j)/cm(j-i);
                cb(j,i)=cb(i,j);
                coeff(i,j)= (N-(j-i)+1)^2;
                coeff(j,i)=coeff(i,j);
            end
        end
        coeff(coeff==0)=nan;
        
        s=[];
        for i=3:N-min([N-2, 40])
            for j=i+1:min([N-2, 40])
                c=cb(i:i+1,j:j+1);
                a=find(c);c=c(a);a=isfinite(c);a=find(a);
                if (size(a,1)==4)
                    s(end+1)=var(c(a))/mean(c(a))^2;
                end
            end
        end
        nrm=median(s); % normalisation factor
        
        % Define distance card for the hierarchical clustering
        card = 1-corr(card0);
        
        % Delete distances computed with too few componant
        nb_min=10;  % min number of values that 2 columns must have in commun
        whereValues = double(card0>0); % no 0 values
        card(whereValues'*whereValues<nb_min)=nan;
        
        % Initialisation of the variables
        dist = diag(card,1); % dist between successive clusters. nan if clusters cannot be merged
        allDomainsLabels=[(1:N)' (1:N)'];
        currentDomainsLabels=(1:N)';
        spreadBetween = diag(card,1);
        domain = [(1:N)' (1:N)'];
        step = 0; % incremented if we remove the boundary
        if ir==1
        Z = zeros(1,3); % matrix to fill. Format = the same that the result of linkage matlab function
        end

        % While dist has no nan value(s) we try to merge domains = to delete boundaries
        while sum(isnan(dist))<length(dist)
            
            % Minimum distance to find the two closest domains to merge
            [m, idxMinDist] = min(dist);
            
            % List of the bins of the 2 domains to merge
            liste1=allDomainsLabels(currentDomainsLabels(idxMinDist),2:end);
            liste1(liste1==0)=[];
            liste2=allDomainsLabels(currentDomainsLabels(idxMinDist+1),2:end);
            liste2(liste2==0)=[];
            liste12 = [liste1 liste2]; liste12 = sort(liste12);
            
            % dt for the new cluster we want to form
            insideCard = card0(liste12, liste12);
            nis = size(insideCard,1);
            cm=zeros(nis-1,1);
            for i=1:nis-1
                d=diag(insideCard,i);a=find(d);d=d(a);a=isfinite(d);a=find(a);
                cm(i)=mean(d(a));
            end
            cb=zeros(nis);
            for i=1:nis
                for j=i+1:nis
                    cb(i,j)=insideCard(i,j)/cm(j-i);
                    cb(j,i)=cb(i,j);
                end
            end
            cbs=cb;
            for i=1:size(cbs,1)
                cbs(i,i)=nan;
            end
            ll1=length(liste1);
            cbs1=cbs(1:ll1,1:ll1); cbs2=cbs(ll1+1:end,ll1+1:end); cbsinter=cbs(ll1+1:end,1:ll1);
            a=find(cb);cb=cb(a);a=isfinite(cb);a=find(a);
            current_dt = var(cb(a))/nrm;
            
            if sigmaZero==0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Conditions to keep or not the boundary - Default case
                if current_dt >= sigmaPlus
                    dist(idxMinDist)=nan;
                    continue;  % We keep this boundary & look for a next boundary to remove
                elseif current_dt <= sigmaMinus
                elseif ( trimmean([cbs1(:); cbs2(:)],25)< trimmean(cbsinter(:),25) )
                elseif ( length(liste1)>=4 && length(liste2)>=4 && length(liste12)>=10) % di test
                    a_di=nan(length(insideCard),1);
                    b_di=nan(length(insideCard),1);
                    if mod(nis,2)==0
                        tmp = [0:nis/2-1 nis/2-1:-1:0];
                    else
                        tmp = [0:(nis-1)/2 (nis-1)/2-1:-1:0];
                    end
                    debut=max(ceil(length(liste1)/2)+1, ll1-10);
                    fin=min(ll1+length(liste2)-ceil(length(liste2)/2), ll1+11);
                    for i=debut:fin
                        a_di(i) = nansum(insideCard(i,i-tmp(i):i-2)); % upstream contacts
                        b_di(i) = nansum(insideCard(i,i+2:i+tmp(i))); % downstream contacts
                    end
                    e_di = nanmean([a_di b_di],2);
                    di = sign(b_di-a_di).*((a_di-e_di).^2./e_di+(b_di-e_di).^2./e_di);
                    diff_rel = abs(a_di-b_di)./e_di;
                    if ( di(ll1)<di(ll1+1) && quantile(1.0*(di(debut+1:ll1)<=0),0.3333)>0 && ...
                            quantile(1.0*(di(ll1+1:fin-1)>=0),0.3333)>0 && ...
                            quantile(1.0*(diff_rel(debut+1:fin-1)>0.1),0.3333)>0)
                        dist(idxMinDist)=nan;
                        continue; % We keep this boundary & look for a next boundary to remove
                    end
                elseif ( length(liste12)>=10 && ( length(liste1)<=4 || length(liste2)<=4 )) % Chaining
                    if length(liste1)>=length(liste2)
                        st=1; en=ll1;
                        adj=log2([insideCard(st:en-1,ll1) insideCard(st+1:en,ll1+1)]);
                        contactB=(adj(:,1)-adj(:,2))./nanmean(adj,2);
                    else
                        st=ll1+1; en=ll1+length(liste2);
                        adj=log2([insideCard(st:en-1,ll1) insideCard(st+1:en,ll1+1)]);
                        contactB=(-adj(:,1)+adj(:,2))./nanmean(adj,2);
                    end
                    
                    if ( quantile((contactB>0)*1.0,0.20)>0 && quantile((contactB>0.01)*1.0,0.50)>0 ) % boundary between a small and a big domain
                        dist(idxMinDist)=nan;
                        continue;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                if current_dt >= sigmaZero(ith)
                    dist(idxMinDist)=nan;
                    continue;  % We keep this boundary & look for a next boundary to remove
                end
            end
            
            step=step+1; % step is incremented if the tests are passed
            for ii=1:nis
                insideCard(ii,ii)=nan;
            end
            
            % Compute the distance between the new cluster and its neighbors
            % If the new cluster has 2 neighbors
            if idxMinDist>1 && idxMinDist<size(currentDomainsLabels,1)-1
                listePrevious=allDomainsLabels(currentDomainsLabels(idxMinDist-1),2:end);
                listePrevious(listePrevious==0)=[];
                listeNext=allDomainsLabels(currentDomainsLabels(idxMinDist+2),2:end);
                listeNext(listeNext==0)=[];
                
                insideCardPrevious = card(listePrevious,listePrevious);
                nPrevious=size(insideCardPrevious,1);
                insideCardNext = card(listeNext,listeNext);
                nNext=size(insideCardNext,1);
                
                nD = max([nis nPrevious nNext]);
                diago = nan(nD,1);
                for ii=1:nis
                    diago(ii)=nanmean(diag(insideCard,ii));
                end
                diagoPrevious = nan(nD,1);
                if nPrevious>1
                    for ii=1:nPrevious
                        diagoPrevious(ii)=nanmean(diag(insideCardPrevious,ii));
                    end
                end
                diagoNext = nan(nD,1);
                if nNext>1
                    for ii=1:nNext
                        diagoNext(ii)=nanmean(diag(insideCardNext,ii));
                    end
                end
                
                betweenCard1 = card(listePrevious,liste12);
                betweenTmpCard1 = nan(nPrevious+nis);
                betweenTmpCard1(1:nPrevious,nPrevious+1:end)=betweenCard1;
                diagoBetween1 = nan(nD,1);
                for ii=1:max([nis nPrevious]) - 1
                    diagoBetween1(ii)=nanmean(diag(betweenTmpCard1,ii));
                end
                coeff1 = coeff(listePrevious,liste12);
                
                betweenCard2 = card(liste12, listeNext);
                betweenTmpCard2 = nan(nis+nNext);
                betweenTmpCard2(1:nis,nis+1:end)=betweenCard2;
                diagoBetween2 = nan(nD,1);
                for ii=1:max([nis nNext]) - 1
                    diagoBetween2(ii)=nanmean(diag(betweenTmpCard2,ii));
                end
                coeff2 = coeff(liste12, listeNext);
                
                spreadBetween(idxMinDist-1) = nansum(nansum(coeff1.*betweenCard1))/nansum(coeff1(:));
                spreadBetween(idxMinDist+1) = nansum(nansum(coeff2.*betweenCard2))/nansum(coeff2(:));
                
                if ~isnan(dist(idxMinDist-1))
                    dist(idxMinDist-1)=spreadBetween(idxMinDist-1);
                end
                
                if ~isnan(dist(idxMinDist+1))
                    dist(idxMinDist+1)=spreadBetween(idxMinDist+1);
                end
                
                % If the new cluster has 1 neighbor
            elseif idxMinDist==1 && idxMinDist~=size(currentDomainsLabels,1)-1
                
                listeNext=allDomainsLabels(currentDomainsLabels(idxMinDist+2),2:end);
                listeNext(listeNext==0)=[];
                insideCardNext = card(listeNext,listeNext);
                nNext=size(insideCardNext,1);
                nD = max([nis nNext]); diago = nan(nD,1);
                for ii=1:nis
                    diago(ii)=nanmean(diag(insideCard,ii));
                end
                diagoNext = nan(nD,1);
                if nNext>1
                    for ii=1:nNext
                        diagoNext(ii)=nanmean(diag(insideCardNext,ii));
                    end
                end
                betweenCard2 = card(liste12, listeNext);
                betweenTmpCard2 = nan(nis+nNext);
                betweenTmpCard2(1:nis,nis+1:end)=betweenCard2;
                diagoBetween2 = nan(nD,1);
                for ii=1:max([nis nNext]) - 1
                    diagoBetween2(ii)=nanmean(diag(betweenTmpCard2,ii));
                end
                coeff2 = coeff(liste12, listeNext);
                
                spreadBetween(idxMinDist+1) = nansum(nansum(coeff2.*betweenCard2))/nansum(coeff2(:));
                if ~isnan(dist(idxMinDist+1))
                    dist(idxMinDist+1)=spreadBetween(idxMinDist+1);
                end
                
                
                % If the new cluster has 1 neighbor
            elseif idxMinDist~=1 && idxMinDist==size(currentDomainsLabels,1)-1
                listePrevious=allDomainsLabels(currentDomainsLabels(idxMinDist-1),2:end);
                listePrevious(listePrevious==0)=[];
                insideCardPrevious = card(listePrevious,listePrevious);
                nPrevious=size(insideCardPrevious,1);
                nD = max([nis nPrevious]);
                diago = nan(nD,1);
                for ii=1:nis
                    diago(ii)=nanmean(diag(insideCard,ii));
                end
                diagoPrevious = nan(nD,1);
                if nPrevious>1
                    for ii=1:nPrevious
                        diagoPrevious(ii)=nanmean(diag(insideCardPrevious,ii));
                    end
                end
                betweenCard1 = card(listePrevious,liste12);
                betweenTmpCard1 = nan(nPrevious+nis);
                betweenTmpCard1(1:nPrevious,nPrevious+1:end)=betweenCard1;
                diagoBetween1 = nan(nD,1);
                for ii=1:max([nis nPrevious]) - 1
                    diagoBetween1(ii)=nanmean(diag(betweenTmpCard1,ii));
                end
                coeff1 = coeff(listePrevious,liste12);
                spreadBetween(idxMinDist-1) = nansum(nansum(coeff1.*betweenCard1))/nansum(coeff1(:));
                if ~isnan(dist(idxMinDist-1))
                    dist(idxMinDist-1)=spreadBetween(idxMinDist-1);
                end
            end
            
            % Update variables after removing the current boundary  
            if ir==1
            Z(step,1:3)=[sort([currentDomainsLabels(idxMinDist), currentDomainsLabels(idxMinDist+1)]) m];
            end
            currentDomainsLabels(idxMinDist)=max(allDomainsLabels(:,1))+1;
            currentDomainsLabels(idxMinDist+1)=[];
            allDomainsLabels(end+1,1:length(liste12)+1)=[max(allDomainsLabels(:,1))+1 liste12];
            spreadBetween(idxMinDist)=[];
            dist(idxMinDist)=[];
            domain(idxMinDist,2)=domain(idxMinDist+1,2);
            domain(idxMinDist+1,:)=[];
        end
        
        % Change value in domain if some columns/rows have been excluded from the
        % study at the begining

        trueIdx=1:Ninitial;
        trueIdx(notEnoughSignal)=[];
        domain(:,1) = trueIdx(domain(:,1));
        domain(:,2) = trueIdx(domain(:,2));
  
        compteur=1; seq_domain=zeros(1,Ninitial);
        for ii=1:size(domain,1)
            seq_domain(domain(ii,1):domain(ii,2))=compteur;
            compteur=compteur+1;
        end
        
        k = [true;diff(notEnoughSignal(:))~=1 ];
        s = cumsum(k);
        nes_dom_size =  histc(s,1:s(end));
        nes_dom_size = nes_dom_size(s);
        
        seq_domain(notEnoughSignal(nes_dom_size>2))=0; % if in one tad, it miss 1 or 2 columns we do not divide the tad in two
        % We delete bins without enought signal if there are 3 or more
        % consecutive
        if seq_domain(end)==0
            if ~isempty(notEnoughSignal) && notEnoughSignal(1)==1
                domain = [[ find(diff(seq_domain) & seq_domain(2:end)~=0)'+1] [find(diff(seq_domain) & seq_domain(1:end-1)~=0)']];
            else
                domain = [[ 1; find(diff(seq_domain) & seq_domain(2:end)~=0)'+1] [find(diff(seq_domain) & seq_domain(1:end-1)~=0)']];
            end
        else
            if ~isempty(notEnoughSignal) && notEnoughSignal(1)==1
                domain = [[ find(diff(seq_domain) & seq_domain(2:end)~=0)'+1] [find(diff(seq_domain) & seq_domain(1:end-1)~=0)'; length(seq_domain)]];
            else
                domain = [[ 1; find(diff(seq_domain) & seq_domain(2:end)~=0)'+1] [find(diff(seq_domain) & seq_domain(1:end-1)~=0)'; length(seq_domain)]];
            end
        end
        
        if ir==1
            % Change value in Z if some columns/rows have been excluded
            % from the study at the begining
            trueIdx=1:2*Ninitial;
            trueIdx(notEnoughSignal)=[];
            Z(:,1)=trueIdx(Z(:,1));
            Z(:,2)=trueIdx(Z(:,2));
            % Fill Z to get the good format for dendrogram
            missingvalues = setdiff(1:2*(Ninitial-1),unique(Z(:,1:2)));
            missingvalues = reshape(missingvalues,[2, length(missingvalues)/2])';
            missingvalues(:,3)=nan;
            Z=[Z; missingvalues];
            
            if sigmaZero==0
                if ( sum(strcmp(def,'sigmaMinus'))==1 && sum(strcmp(def,'sigmaPlus'))==1 )
                    nameOptionHierarchy='';
                else
                    nameOptionHierarchy=['_sigmaMinus' num2str(sigmaMinus) '_sigmaPlus' num2str(sigmaPlus)];
                end
                if saveDomains==1
                    dlmwrite([Path name nameOptionHierarchy '_domains.txt'],domain,'delimiter', '\t','precision',8);
                end
                if saveZ==1
                    dlmwrite([Path name nameOptionHierarchy '_Z.txt'],Z,'delimiter', '\t','precision',8);
                end
            else
                nameOptionHierarchy = ['_sigmaZero' num2str(sigmaZero(ith))];
                if saveDomains==1
                    dlmwrite([Path name nameOptionHierarchy '_domains.txt'],domain,'delimiter', '\t','precision',8);
                end
                if saveZ==1
                    dlmwrite([Path name nameOptionHierarchy 'Z.txt'],Z,'delimiter', '\t','precision',8);
                end
            end
                      
            
            dom = domain; % domain can be changed after if nr>1
            
            Zoutput(:,:,ith) = Z;
            domoutput(1:size(dom,1),1:2,ith) = dom;
            
            
        end
        
        % update pd & pb
        for ii=1:size(domain,1)
            pd(domain(ii,1):domain(ii,2),domain(ii,1):domain(ii,2)) = pd(domain(ii,1):domain(ii,2),domain(ii,1):domain(ii,2))+1;
        end
        b = [b; unique([domain(:,1)-0.5; domain(:,2)+0.5])];
        
    end % end resampling
    [y,x] = hist(b,0.5:Ninitial+0.5);
    pb=[x(:), y(:)/nr];
    pd = pd./nr;
    pdoutput(:,:,ith) = pd;
    pboutput(1:size(pb,1),1:2,ith) = pb;
    
    if nr>1
        nameOptionResampling=['_nr' num2str(nr) '_nl' num2str(nl)];
        if savePD==1
            dlmwrite([Path name nameOptionHierarchy nameOptionResampling '_pd.txt'],pd,'delimiter', '\t','precision',8);
        end
        if savePB==1
            dlmwrite([Path name nameOptionHierarchy nameOptionResampling '_pb.txt'],pb,'delimiter', '\t','precision',8);
        end
    else
        nameOptionResampling = '';
    end
    
    if plotFigures==1 || saveFigures==1
        if nr>1
            [fig1, fig2, fig3]=pIC_Finder(Cinitial,dom,'Z',Z,'PD',pd,'PB',pb,'Name',[nameOptionHierarchy nameOptionResampling]);
        else
            [fig1, fig2, fig3]=pIC_Finder(Cinitial,dom,'Z',Z,'PD',0,'PB',[0 0],'Name',[nameOptionHierarchy nameOptionResampling]);
        end
        
        if saveFigures==1
            print(fig1,[Path name nameOptionHierarchy nameOptionResampling '_fig1.pdf'],'-dpdf');
            print(fig2,[Path name nameOptionHierarchy nameOptionResampling '_fig2.pdf'],'-dpdf');
            if ~isempty(fig3)
            print(fig3,[Path name nameOptionHierarchy nameOptionResampling '_fig3.pdf'],'-dpdf');
            end
        end
        if plotFigures==0
            close(fig1); close(fig2); close(fig3);
        end
    end
        
end % end hierarchy (different threshold)

Z = Zoutput;
dom=domoutput;
pd=pdoutput;
pb=pboutput;


function [fig1, fig2, fig3]=pIC_Finder(HiCmap,dom, varargin)
%
%   by Noelle HADDAD, Cédric VAILLANT and Daniel JOST
%   May 2016
%   Please read the LICENSE file before using IC_Finder
%   Any remarks and suggestions should be addressed to Daniel JOST daniel.jost@imag.fr
%
%pIC_Finder(HICMAP,DOM) plots an HiCmap with domain boundaries.
%   HICMAP is a symmetric square matrix or the path to a .txt or .mat file 
%   containing it. DOM is a two-columns matrix or DOM can be the path to a
%   .txt or .mat file containing a two-columns matrix. In this txt
%   files there must be no header. Each line must represent one domain, in
%   the first and second columns, there are respectively the first and last
%   bins of the domain.
%
%   pIC_Finder(HICMAP,DOM, 'PARAM1',val1, 'PARAM2',val2, ...)
%   specifies optional parameter name/value pairs to control plots.
%   Parameters are : 
%
%   'Z'      -   3-columns matrix (or path to a .txt or .mat file 
%                containing it) enabling to plot in a new figure the 
%                dendrogram of the hierarchical clustering. Z is an output 
%                of IC_finder function. If Z is not specified, no 
%                dendrogram plot.
%
%   'PD'     -   Symmetric square matrix (or path to a .txt or .mat file
%                containing it) representing the probability that two bins
%                belongs to the same domain. PD is an output of IC_finder
%                function. If you specify PD parameter, this matrix will be
%                plotted in a new figure.
%
%   'PB'     -   2-columns matrix (or path to a .txt or .mat file 
%                containing it) representing the probability that a bin is 
%                detected as a boundary. PB is an output of IC_finder
%                function. If you specify PB parameter, the probability to
%                be a predicted as a boundary as a function of the bin will
%                be plotted in a new figure or in the one where there is PD
%                if PD has been specified.
%
%   'Name'       - Name to save images. By default, if HICMAP is a string,
%                Name is equal to the HICMAP filename without the extension.
%                If HICMAP is not a path but a matrix Name is by default :
%                image_YearMonthDay_Hour:Minute:Second .
%
%   'START'  -   Float specifying the beginning of the region one wants to
%                plot (in bin). By default, START=1.
%
%   'EN'     -   Float specifying the end of the region one wants to plot
%                (in bin). By default , END=size of HICMAP.
%
%
%   Examples :
%
%      % To plot everything, ie. the HiC map, domains, dendrogram, PD and PB
%      f='example/maptest_1'; %filename
%      pIC_Finder([f '.txt'], [f '_domains.txt'], 'Z', [f '_Z.txt'],...
%      'PD',[f '_pd_nr100_nl1.txt'],'PB',[f '_pb_nr100_nl1.txt'],'start',1,'end',80);
%
%      % To plot the HiC map, domains and PD :
%      pIC_Finder([f '.txt'], [f '_domains.txt'],'PD',[f '_pd_nr100_nl1.txt']);
%
%   See also IC_FINDER
%


% Inputs
if ischar(HiCmap)
    card = load(HiCmap); 
    [pathstr,defaultName] = fileparts(HiCmap);
    if ~isempty(pathstr)
        defaultPath=[pathstr '/'];
    else
        defaultPath=pathstr;
    end
else
    card = HiCmap;
    defaultPath='';
    current_date = clock;
    name=['card_' num2str(current_date(1)) num2str(current_date(2),'%02.0f') num2str(current_date(3),'%02.0f')...
        '_' num2str(current_date(4),'%02.0f') ':' num2str(current_date(5),'%02.0f') ':' num2str(current_date(6),'%02.0f')];
    defaultName=strrep(name,' ','_');
    nameD = 1; % if name is date
end

N=size(card,1);
if ischar(dom)
    domain = load(dom);
else
    domain = dom;
end


p = inputParser;
p.FunctionName = 'pIC_Finder';

defaultStart = 1;
defaultEn = N;
defaultZ=0; 
defaultPD=0; 
defaultPB=[0 0]; 

addRequired(p,'HiCmap');
addRequired(p,'dom');
addParamValue(p,'Z',defaultZ)
addParamValue(p,'PD',defaultPD)
addParamValue(p,'PB',defaultPB)
addParamValue(p,'name',defaultName)
addParamValue(p,'start',defaultStart,@isnumeric)
addParamValue(p,'en',defaultEn,@isnumeric)

parse(p,HiCmap,dom,varargin{:})

def = p.UsingDefaults;

st = p.Results.start;
en = p.Results.en;
Z = p.Results.Z;
pd = p.Results.PD;
pb = p.Results.PB;
name = p.Results.name;

if ~sum(strcmp(def,'Z'))==1
    if ischar(Z)
        Z = load(Z);
    else
        validateattributes(Z,{'numeric'},{'ncols',3},'pIC_Finder','Z')
    end
end
if ~sum(strcmp(def,'PD'))==1
    if ischar(pd)
        pd = load(pd);
    else
        validateattributes(pd,{'numeric'},{'square'},'pIC_Finder','PD')
    end
end
if ~sum(strcmp(def,'PB'))==1
    if ischar(pb)
        pb = load(pb);
    else
        validateattributes(pb,{'numeric'},{'ncols',2},'pIC_Finder','PB')
    end
end

validateattributes(st,{'numeric'},{'scalar','>=',0,'<=',N})
validateattributes(en,{'numeric'},{'scalar','>=',st,'<=',N})
validateattributes(name,{'char'},{})

% 1 full segmentation
fig1=figure; set(gcf,'color','w','units','normalized','outerPosition',[0 0 1 1],'Name',name)
colormap(jet)
subplot('position', [0.2942 0.11 0.4466 0.8150])
imagesc(log2(card)); hold on; 
for i=1:length(domain)
    plot([1 N], [domain(i,1)-0.5 domain(i,1)-0.5], 'w-','LineWidth',2)
    plot([1 N], [domain(i,2)+0.5 domain(i,2)+0.5], 'w-','LineWidth',2)
    plot([domain(i,1)-0.5 domain(i,1)-0.5], [1 N], 'w-','LineWidth',2)
    plot([domain(i,2)+0.5 domain(i,2)+0.5], [1 N], 'w-','LineWidth',2)
end
axis([st en st en])
title('HiC map with domains')

% 1 semi segmentation & the dendrogram
fig2=figure; set(gcf,'color','w','units','normalized','outerPosition',[0 0 1 1],'Name',name);
colormap(jet)
if ~isequal(Z,0)
    subplot('position', [0.13 0.2121 0.3347 0.6108]);
else
    subplot('position', [0.2942 0.11 0.4466 0.8150]);
end
imagesc(log2(card)); hold on; 
ax=gca;
xv = [0.5 N+0.5 N+0.5];
yv = [0.5 0.5 N+0.5]; 
for i=1:length(domain)
    xplot = 1:N;
    yplot = 0*xplot+domain(i,1)-0.5;
    yplot2 = 0*xplot+domain(i,2)+0.5;
    in = inpolygon(xplot,yplot,xv,yv);
    in2 = inpolygon(xplot,yplot2,xv,yv);
    plot(xplot(in), yplot(in), 'w-','LineWidth',2)
    plot(xplot(in2), yplot2(in2), 'w-','LineWidth',2)
    
    yplot = 1:N;
    xplot = 0*yplot+domain(i,1)-0.5;
    xplot2 = 0*yplot+domain(i,2)+0.5;
    in = inpolygon(xplot,yplot,xv,yv);
    in2 = inpolygon(xplot2,yplot,xv,yv);
    plot(xplot(in), yplot(in), 'w-','LineWidth',2)
    plot(xplot2(in2), yplot(in2), 'w-','LineWidth',2)
end
axis([st en st en])
title('HiC map with domains')

if ~isequal(Z,0)
    subplot('position',[0.49 0.2121 0.3347 0.6108]); 
    %dendrogram(Z,N,'Reorder',1:N,'Orientation','right');
    
    %%%%% dendrogram(Z,N,'Reorder',1:N,'Orientation','right'); %%%%%
    % plot manually because dendrogram function not available with octave
    m = size (Z,1);
    n = m + 1;
    t = (1:m)'; 
    nc = max(max(Z(:,1:2)));
    
    p = zeros (nc,2); x = zeros (m,2);
    p(1:nc,1)=1:nc;
    for i = 1:m
        p(n+i,1)   = mean (p(Z(i,1:2),1));
        p(n+i,2)   = Z(i,3);
        x(i,1:2) = p(Z(i,1:2),1);
    end
    
    line (Z(:,[3 3])',x','color','k');
    
    [~,tf]  = ismember (1:nc, Z(:,1:2));
    [ind,~] = ind2sub (size (Z(:,1:2)), tf);
    y       = [p(1:nc,2) Z(ind,3)];
    line (y',[p(1:nc,1) p(1:nc,1)]','color','k');
    
    yticks = 1:n;
    yl_txt = arrayfun (@num2str, (1:n)','uniformoutput',false);
    set (gca,'yticklabel',yl_txt,'ytick',yticks);
    axis ([0 max(Z(:,3))+0.1*min(Z(:,3)) 0.5 n+0.5]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ylim([st en])
    set(gca,'YDir','reverse')
    ax(2)=gca; set(ax(2),'yticklabel',[]);
    linkaxes(ax,'y')
    title('Hierarchical clustering dendrogram')
end


% Segmentation and statistical outputs
fig3=[];
if ( ~isequal(pd,0) || ~isequal(pb,[0 0]) )
    fig3=figure; set(gcf,'color','w','units','normalized','outerPosition',[0 0 1 1],'Name',name)
    colormap(jet);
    if ( ~isequal(pd,0) )
        subplot('position',[0.13 0.2121 0.3347 0.6108])
    else
        subplot('position', [0.2942 0.11 0.4466 0.8150])
    end
    imagesc(log2(card));  hold on; ax=gca;
    set(ax,'xticklabel',[]);
    for i=1:length(domain)
        plot([1 N], [domain(i,1)-0.5 domain(i,1)-0.5], 'w-')
        plot([1 N], [domain(i,2)+0.5 domain(i,2)+0.5], 'w-')
        plot([domain(i,1)-0.5 domain(i,1)-0.5], [1 N], 'w-')
        plot([domain(i,2)+0.5 domain(i,2)+0.5], [1 N], 'w-')
    end
    axis([st en st en])
    title('HiC map and domains')
    
    if ( ~isequal(pb,[0 0]) )
    subplot('position', [0.13 0.1461 0.3347 0.0580]); 
    plot(pb(:,1),pb(:,2)); ax(end+1)=gca;
    linkaxes(ax(1:2),'x')
    axis([st en 0 1])
    ylabel('P_B')
    end
    
    if ( ~isequal(pd,0) )
    subplot('position',[0.51 0.2121 0.3347 0.6108]); imagesc(pd); colormap(jet)
    ax(end+1)=gca;
    setappdata(ax(1), 'XLim_listeners', linkprop([ax(1),ax(end)],'XLim'));
    setappdata(ax(1), 'YLim_listeners', linkprop([ax(1),ax(end)],'YLim'));
    title('P_D')
    end
end
