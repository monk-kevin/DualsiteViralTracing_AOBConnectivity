function fcnProbabilisticToyModels(nLbl,nBoots)
%% Description
% This function creates a visualization for the estimation analysis with
% probabilistic, ground-truth models. Namely, it visualizes the ranges of
% co-labeling rates one could observe for a given ground-truth
% co-projection rate. This relationship will vary based on the number of
% labeled neurons observed and this function shows that relationship by
% varying the observed sample size.
%% Inputs
% nLbl: the sample size of a simulated dataset corresponding to the number
% of total labeled neurons. This parameter will be used to show the range
% of co-projection probabilities an nLbl-sized population falls within.
% Default value is 100.
%
% nBoots: the number of iterations performed when establishing the
% ground-truth, probabilistic models. Default value is 1000.
%% Output
% A figure demonstrating the range of co-projection probabilities a sample
% of a given size is consistent with. The top row shows extended results
% for one specific sample size (defined by nLbl) and bottom shows range
% estimates across simulated datasets of varying sample sizes.

%% ground truth toy models
% set up for probabilistic models
if nargin == 0
    nLbl = 100; % default sample size for top row illustration
    nBoots = 1000; %default bootstraps for model creation
elseif nargin == 1
    nBoots = 1000;
end
GroundTruthProp = 0:0.01:1;
ChancePropLbl = nan(nBoots,numel(GroundTruthProp));

%figure set up
figure
FS = 20; MS = 16;


for pp = 1:numel(GroundTruthProp)
    for nn = 1:nBoots
        RandLbl = rand(nLbl,1);
        ChancePropLbl(nn,pp) = sum(RandLbl<=GroundTruthProp(pp))/nLbl;
    end
end

%visualize data from n = 100 with exemplar proportions to show range of
%consistent co-projection rates
ExampleProps = [25 50 75];

CLR = {'r','g','b'};
subplot(211);hold on
for ii = 1:numel(ExampleProps)
    histogram(ChancePropLbl(:,ExampleProps(ii)+1),0:0.01:1,'normalization','probability','facecolor',CLR{ii},'facealpha',1)
    xline(prctile(ChancePropLbl(:,ExampleProps(ii)+1),[5 95]),'color',CLR{ii},'linestyle','--','linewidth',2)
end
set(gca,'fontsize',FS)
xlabel('Observed Co-Projection Probability');ylabel('Proportion of Iterations');

% varying neural population size to show ranges for varying sample size
PopSize = [5 25 50 100 500 1000];
for ii = 1:numel(PopSize)
    nLbl = PopSize(ii);
    for pp = 1:numel(GroundTruthProp)
        for nn = 1:nBoots
            RandLbl = rand(nLbl,1);
            ChancePropLbl(nn,pp) = sum(RandLbl<=GroundTruthProp(pp))/nLbl;
        end
    end
    PRTIL = [prctile(ChancePropLbl,5);prctile(ChancePropLbl,95)];
    subplot(2,numel(PopSize),ii + numel(PopSize))
    fill([GroundTruthProp fliplr(GroundTruthProp)],...
        [PRTIL(1,:) fliplr(PRTIL(2,:))],...
        ones(1,3)*.6)
    line([0 1],[0 1],'color','k','linestyle','--','linewidth',2)
    title(['n = ' num2str(PopSize(ii))])
    axis([0 1 0 1])
    xlabel({'Ground Truth';'Co-Projection Proportion'})
    ylabel('Measured Co-Projection Proportion')
    set(gca,'fontsize',FS)
end

