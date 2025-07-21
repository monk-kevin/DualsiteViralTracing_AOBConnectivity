function fcnPlotPropProj_CrossFluorophore(AniSliceTbl_avg)
%% Description
% This dataset contains neurons that were retrogradely labeled using a
% virus that contained either GFP or tdTomato. This function plots the
% proportion of MCs labeled with either flurophore for a given area to
% determine whether the flurophore used creates a difference in labeling
% probabilities
%% Input
% AniSliceTbl_avg: triplicate-averaged data for each experimental unit
% (animal/hemisphere pair). Created using the fcnCreateAniSliceTbl function
%% Output
% figure of paired bar graphs demonstrating the proportion of AOB MCs that
% are labeled as projecting to a given area split by the fluroophore used
% in that specific labeling experiment.
%%
% divide table into dual and single site cohorts
AREAs = unique([AniSliceTbl_avg.GFPSource;AniSliceTbl_avg.tdTomSource]);

propGFP = cell(1,numel(AREAs));
propTdTom = cell(1,numel(AREAs));
MS = 24;

for ii = 1:numel(AREAs)
    propGFP{ii} = AniSliceTbl_avg.propGFP(strcmp(AniSliceTbl_avg.GFPSource,AREAs{ii}));
    propTdTom{ii} = AniSliceTbl_avg.propTdTom(strcmp(AniSliceTbl_avg.tdTomSource,AREAs{ii}));
end

figure
hold on
b = bar([cellfun(@median,propGFP);cellfun(@median,propTdTom)]');
b(1).FaceColor = [0 1 0]*.7;
b(2).FaceColor = [1 0 0]*.7;

for ii = 1:numel(AREAs)
    plot(ones(numel(propGFP{ii}),1)*b(1).XEndPoints(ii),propGFP{ii},'k.','markersize',MS)
    plot(ones(numel(propTdTom{ii}),1)*b(2).XEndPoints(ii),propTdTom{ii},'k.','markersize',MS)

    % defining the proportion of labeled neurons based on fluorophore used
    X = b(1).XEndPoints;Y1 = cellfun(@median,propGFP);ERR = cellfun(@std,propGFP);N = cellfun(@numel,propGFP);
    errorbar(X,Y1,ERR./sqrt(N),'color','k','linestyle','none','linewidth',4)
    
    X = b(2).XEndPoints;Y2 = cellfun(@median,propTdTom);ERR = cellfun(@std,propTdTom);N = cellfun(@numel,propTdTom);
    errorbar(X,Y2,ERR./sqrt(N),'color','k','linestyle','none','linewidth',4)
    
    % stats
    p(ii) = ranksum(propGFP{ii},propTdTom{ii});
    line([b(1).XEndPoints(ii) b(2).XEndPoints(ii)],[0.575 0.575],'color','k','linewidth',2)
    text(b(1).XEndPoints(ii),0.6,['p = ' num2str(round(p(ii),2))],'fontsize',18)
    
end
xlabel('Target Area');
ylabel('Proportion of Back-Labeled MCs');
set(gca,'xtick',1:3,'xticklabel',AREAs,'fontsize',16)
axis([.6 3.4 0 0.65])