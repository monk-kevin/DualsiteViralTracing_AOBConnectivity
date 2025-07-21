function [AniSliceTbl,AniSliceTbl_avg] = fcnCreateAniSliceTbl(XLS_IDinfo,path2Results)
%% Description
% This function returns two tables that contains all relevant information
% for connectivity analysis based on the animal information provided and
% the appropriate paths to saved cell counter results
%% Inputs
% XLS_IDinfo: path to a master spreadsheet with animal metadata
% path2Results: path to cell counter results created in ImageJ that
% contains the number and spatial location of each labeled neuron
%% Outputs
% AniSliceTbl: a table that contains animal identifying information, virus
% strategy, the number of labeled neurons for each label, as well as the
% spatial organization within the AOB
%
% AniSliceTbl_avg: a table that averages across triplicate sections within
% AniSliceTbl variable. 
%%
% read master info table to get complete list of injected animals
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
tbl_ID = readtable(XLS_IDinfo);

%define path to cell results and create mask to look for results files
strMask = '\*.csv';
for ii = 1:numel(path2Results)
    tempD = dir([path2Results{ii} strMask]);
    if ii == 1
        D = tempD;
    else
        D = [D;tempD];
    end
end

% Cell Counter Categories:
%1: neurotrace/blue
%2: GFP/green
%3: tdTomato/red
%4: GFP + tdTomato OVERLAP
%5: ANTERIOR LANDMARK
%6: POSTERIOR LANDMARK

AniSliceTbl = table;
for ii = 1:numel(D)
    if rem(ii,5)==0
        disp(['Starting Slice #' num2str(ii) ' of ' num2str(numel(D))])
    end
    if length(D(ii).name) == 31
        tAniID = D(ii).name(13:18);
        if strcmp(tAniID,'39098N')
            AniID = 390981;
        elseif strcmp(tAniID,'39098L')
            AniID = 390982;
        elseif strcmp(tAniID,'39098R')
            AniID = 390983;
        elseif strcmp(tAniID,'39099N')
            AniID = 390991;
        elseif strcmp(tAniID,'39100N')
            AniID = 391001;
        elseif strcmp(tAniID,'39100L')
            AniID = 391002;
        elseif strcmp(tAniID,'39100R')
            AniID = 391003;
        elseif strcmp(tAniID,'39101N')
            AniID = 391011;
        elseif strcmp(tAniID,'39101L')
            AniID = 391012;
        elseif strcmp(tAniID,'39101R')
            AniID = 391013;
        elseif strcmp(tAniID,'38123N')
            AniID = 381231;
        end
        %aniID 122551 is 122551L
    elseif length(D(ii).name) == 32 || length(D(ii).name) == 33
        AniID = str2double(D(ii).name(13:18));
    elseif length(D(ii).name) == 34 || length(D(ii).name) == 36
        AniID = str2double(D(ii).name(20:22));
    elseif length(D(ii).name) == 30
        AniID = str2double(D(ii).name(13:16));
    else
        AniID = str2double(D(ii).name(13:15));
    end
    if contains(D(ii).name,'R AOB')||contains(D(ii).name,'RAOB')
        HemiID = 'R';
    else
        HemiID = 'L';
    end
    SliceNum = str2double(D(ii).name(end-4));
    AniSex = tbl_ID.Sex{tbl_ID.AnimalID==AniID};
    GFPSource = tbl_ID.GFPSource{tbl_ID.AnimalID==AniID};
    tdTomSource = tbl_ID.tdTomatoSource{tbl_ID.AnimalID==AniID};
    t = readtable([D(ii).folder '\' D(ii).name]);
    CellCounts=histcounts(t.Type,0.5:4.5);

    tMinDV = t(t.Y==min(t.Y(t.Type==1)),:);
    tMaxDV = t(t.Y==max(t.Y(t.Type==1)),:);
    DV = [tMinDV.X(1) tMinDV.Y(1);tMaxDV.X(1) tMaxDV.Y(1)];

    tNT = t(t.Type==1|t.Type==5|t.Type==6,:);
    tG = t(t.Type==2|t.Type==5|t.Type==6,:);
    tR = t(t.Type==3|t.Type==5|t.Type==6,:);
    tO = t(t.Type==4|t.Type==5|t.Type==6,:);

    %         1st column 0 is more POSTERIOR and 1 is more ANTERIOR.
    %         2nd column 0 is more DORSAL and 1 is more VENTRAL
    [FracAntPost_NT,RotMatStruct] = fcnFindFracAntPost_v2(tNT);
    FracAntPost_g = fcnFindFracAntPost_v2(tG,RotMatStruct);
    FracAntPost_r = fcnFindFracAntPost_v2(tR,RotMatStruct);
    FracAntPost_o = fcnFindFracAntPost_v2(tO,RotMatStruct);


    tempCell = {'AniID','AniSex','HemiID','SliceNum','GFPSource','tdTomSource',...
        'nNT','nGFP','nTdTom','nOverlap','propGFP','propTdTom','propOverlap'...
        'FracAntPostNT','FracAntPostGFP','FracAntPostTdTom','FracAntPostOverlap';...
        AniID,AniSex,HemiID,SliceNum,GFPSource,tdTomSource,...
        CellCounts(1),CellCounts(2),CellCounts(3),CellCounts(4),...
        CellCounts(2)/CellCounts(1),CellCounts(3)/CellCounts(1),CellCounts(4)/CellCounts(1)...
        FracAntPost_NT,FracAntPost_g,FracAntPost_r,FracAntPost_o};
    tempTable = cell2table(tempCell(2,:));tempTable.Properties.VariableNames = tempCell(1,:);
    AniSliceTbl = [AniSliceTbl;tempTable];
end
clear temp*
% averaging across triplicate sections
% note proportion are averages of proportions from above, may wish to calc
% later as avg counts/avg nt (roughly equal numbers, but slightly diff)
ANIs = unique(AniSliceTbl.AniID);
AniSliceTbl_avg = table;
nSteps = 100;
HCbin = 0:1/nSteps:1;
fcnHC_AP = @(X) histcounts(X(:,1),HCbin,'Normalization','probability');
fcnHC_DV = @(X) histcounts(X(:,2),HCbin,'Normalization','probability');
for ii = 1:numel(ANIs)
    t = AniSliceTbl(AniSliceTbl.AniID==ANIs(ii),:);
    HEMIs = unique(t.HemiID);
    for jj = 1:numel(HEMIs)
        t1 = t(strcmp(t.HemiID,HEMIs{jj}),:);
        HCNTap = cellfun(fcnHC_AP,t1.FracAntPostNT,'UniformOutput',false);
        HCgAP = cellfun(fcnHC_AP,t1.FracAntPostGFP,'UniformOutput',false);
        HCtAP = cellfun(fcnHC_AP,t1.FracAntPostTdTom,'UniformOutput',false);
        HCoAP = cellfun(fcnHC_AP,t1.FracAntPostOverlap,'UniformOutput',false);

        HCNTdv = cellfun(fcnHC_DV,t1.FracAntPostNT,'UniformOutput',false);
        HCgDV = cellfun(fcnHC_DV,t1.FracAntPostGFP,'UniformOutput',false);
        HCtDV = cellfun(fcnHC_DV,t1.FracAntPostTdTom,'UniformOutput',false);
        HCoDV = cellfun(fcnHC_DV,t1.FracAntPostOverlap,'UniformOutput',false);
        tempCell = {'AniID','AniSex','HemiID','SliceNum','GFPSource','tdTomSource',...
            'nNT','nGFP','nTdTom','nOverlap','propGFP','propTdTom','propOverlap'...
            'AntPostHC_NT','AntPostHC_GFP','AntPostHC_tdTom','AntPostHC_Overlap'...
            'DorsVentHC_NT','DorsVentHC_GFP','DorsVentHC_tdTom','DorsVentHC_Overlap';...
            ANIs(ii),unique(t1.AniSex),HEMIs{jj},'AVG',unique(t1.GFPSource),unique(t1.tdTomSource),...
            mean(t1.nNT),mean(t1.nGFP),mean(t1.nTdTom),mean(t1.nOverlap),mean(t1.propGFP),mean(t1.propTdTom),mean(t1.propOverlap),...
            mean(cat(1,HCNTap{:}),1,'omitnan'),mean(cat(1,HCgAP{:}),1,'omitnan'),mean(cat(1,HCtAP{:}),1,'omitnan'),mean(cat(1,HCoAP{:}),1,'omitnan'),...
            mean(cat(1,HCNTdv{:}),1,'omitnan'),mean(cat(1,HCgDV{:}),1,'omitnan'),mean(cat(1,HCtDV{:}),1,'omitnan'),mean(cat(1,HCoDV{:}),1,'omitnan')};
        
        tempTable = cell2table(tempCell(2,:));tempTable.Properties.VariableNames = tempCell(1,:);
        AniSliceTbl_avg = [AniSliceTbl_avg;tempTable];
    end
end