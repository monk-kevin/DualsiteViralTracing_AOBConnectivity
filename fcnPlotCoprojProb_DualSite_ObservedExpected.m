function [p,stats] = fcnPlotCoprojProb_DualSite_ObservedExpected(AniSliceTbl_avg)
%% Description
% Plots the observed and expected probabilities of co-labeling for each
% pair of downstream areas
%% Inputs
% AniSliceTbl_avg: triplicate-averaged data for each experimental unit
% (animal/hemisphere pair). Created using the fcnCreateAniSliceTbl function
%% Outputs
% Figure that shows observed (light gray) and expected (dark gray)
% co-labeled probabilities for each pair of area. Expected here is defined
% as the product of the projection probabilities for each area,
% independently.
% p and stats are results of the kruskal-wallis test across distributions
% and values are printed on figure.
%%
AREAs = unique([AniSliceTbl_avg.GFPSource;AniSliceTbl_avg.tdTomSource]);
DualSiteAvg = AniSliceTbl_avg(~strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);
ObsExp = cell(numel(AREAs));
PropOverlap = cell(numel(AREAs));


% for each experimental unit, for each pair of areas, define the expected 
% co-projection proportion as the product of the singular projection 
% probabilities
for ii = 1:numel(AREAs)
    for jj = 1:numel(AREAs)
        if jj ~= ii
            tDual = [DualSiteAvg(strcmp(DualSiteAvg.GFPSource,AREAs{ii})&strcmp(DualSiteAvg.tdTomSource,AREAs{jj}),:);...
                DualSiteAvg(strcmp(DualSiteAvg.tdTomSource,AREAs{ii})&strcmp(DualSiteAvg.GFPSource,AREAs{jj}),:)];
            ObsExp{ii,jj}(:,1) = tDual.propOverlap;
            ObsExp{ii,jj}(:,2) = tDual.propGFP.*tDual.propTdTom;
        end
    end
end
figure
hold on
FS = 16; MS = 24;
ObsExpPairs = {ObsExp{1,2} ObsExp{1,3} ObsExp{2,3}};
for ii = 1:numel(ObsExpPairs)
    b = bar(ii,mean(ObsExpPairs{ii}));b(1).FaceColor = ones(1,3)*.8;b(2).FaceColor = ones(1,3)*.4;
    plot([b(1).XEndPoints b(2).XEndPoints],ObsExpPairs{ii}','color',ones(1,3)*0.5,'linewidth',2)
    errorbar(b(1).XEndPoints,mean(ObsExpPairs{ii}(:,1)),std(ObsExpPairs{ii}(:,1))/sqrt(numel(ObsExpPairs{ii}(:,1))),'color','k','linewidth',3)
    errorbar(b(2).XEndPoints,mean(ObsExpPairs{ii}(:,2)),std(ObsExpPairs{ii}(:,2))/sqrt(numel(ObsExpPairs{ii}(:,1))),'color','k','linewidth',3)
    
    pSR = signrank(ObsExpPairs{ii}(:,1),ObsExpPairs{ii}(:,2));
    line([b(1).XEndPoints b(2).XEndPoints],[0.4 0.4],'color','k','linewidth',3)
    text(b(1).XEndPoints-0.2,0.425,['p = ' num2str(round(pSR,3))],'fontsize',20)
end
axis([0.25 3.75 0 0.5])
TTL = {'BNST+MeA';'BNST+PMCo';'MeA+PMCo'};
set(gca,'xtick',1:3,'xticklabel',TTL,'fontsize',FS)
ylabel('Proportion Co-Labeled MCs')

% stats to be outputted, if required
statsMat = [ObsExpPairs{1}(:,1) ones(size(ObsExpPairs{1}(:,1)));...
    ObsExpPairs{2}(:,1) ones(size(ObsExpPairs{2}(:,1)))*2;...
    ObsExpPairs{3}(:,1) ones(size(ObsExpPairs{3}(:,1)))*3];
[p,~,stats] = kruskalwallis(statsMat(:,1),statsMat(:,2),'off');