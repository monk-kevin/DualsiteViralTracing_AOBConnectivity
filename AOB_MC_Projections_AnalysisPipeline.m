%% Description
% This script will load in cell count data from the Davison Lab server
% and perform analyses to identify connection patterns between the
% accessory olfactory bulb (AOB) and downstream limbic targets (i.e., the
% BNST, MeA, and PMCo). 
%
% Each section calls a function that processes and visualizes data to
% address a specific question. Where noted, the specific function will also
% perform statistical analyses to define significant differences among
% tested distributions.

%% Load in and format raw data into table
XLS_IDinfo = 'Z:\bdibenedictis\BD1 Dual Virus Project\BD InjectedAnimalMasterSpreadsheet_240503_1.xlsx';
path2Results = {...
    'Z:\bdibenedictis\Recounted Data\CellCounterResults_Summer20232';...
    'Z:\bdibenedictis\Recounted Data\CellCountingResults_winter2023_2024';...
    'Z:\bdibenedictis\Recounted Data\CellCountingResults240423'...
    };

[AniSliceTbl,AniSliceTbl_avg] = fcnCreateAniSliceTbl(XLS_IDinfo,path2Results);

%% Projection Probability Analyses to Single Sites
% 1) Plot the proportion of MCs that project to a given area and their
% anatomical distribution within the AOB, collapsed across sex:
SplitSexFlag = 0;
fcnPlotPropProj_AnatDistribution(AniSliceTbl,AniSliceTbl_avg,SplitSexFlag)

% 2) Plot the proportion of MCs that project to a given area based on the
% fluorophore that was injected into that area
fcnPlotPropProj_CrossFluorophore(AniSliceTbl_avg)

% 3) Plot the proportion of MCs that project to a given area split across
% animal sex
SplitSexFlag = 1;
fcnPlotPropProj_AnatDistribution(AniSliceTbl,AniSliceTbl_avg,SplitSexFlag)

%% Co-Projection Probability Analyses to Pairs of Sites
% 1) Probability of a labeled neuron being co-labeled for animals that
% received mixed virus (i.e., both fluorophores) in a single area
fcnPlotCoprojProb_SingleSite(AniSliceTbl_avg)

% 2) Probability of MCs being co-labeled across pairs of areas and compared
% observed probability with the co-labeled probabilities expected if
% neurons were independently projecting to downstream areas.
fcnPlotCoprojProb_DualSite_ObservedExpected(AniSliceTbl_avg); %suppress stats output

%% Defining Significant Co-Labeling Events Through Bootstrap Statistics
nBoots = 1000; %size of the bootstrap population
AlphaPrctile = 99; %statistical threshold defined as Nth percentile
% 1) Dual site cohort:
fcnBootstrapComparison_DualSite(AniSliceTbl_avg,nBoots,AlphaPrctile)

% 2) Single site cohort:
fcnBootstrapComparison_SingleSite(AniSliceTbl_avg,nBoots,AlphaPrctile)

%% Estimating actual co-projection probabilities with probabilistic models
nBoots = 1000; %size of bootstrap population for estimation

% 1) Toy models and visualizing simulated data
fcnProbabilisticToyModels(100, nBoots);

% 2) Dual site cohort:
fcnEstimateTrueProjProb_DualSite(AniSliceTbl_avg, nBoots)

% 3) Single site cohort:
fcnEstimateTrueProjProb_SingleSite(AniSliceTbl_avg, nBoots)
