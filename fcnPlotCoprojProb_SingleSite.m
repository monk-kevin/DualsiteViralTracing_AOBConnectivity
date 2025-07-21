function fcnPlotCoprojProb_SingleSite(AniSliceTbl_avg)
%% Description
% determines the proportion of overlap within the single-site cohort (i.e.,
% animals that received mixed virus in a given downstream area
%% Inputs
% AniSliceTbl_avg: triplicate-averaged data for each experimental unit
% (animal/hemisphere pair). Created using the fcnCreateAniSliceTbl function
%% Outputs
% figure that shows the proportion of labeled MCs that are co-labeled with
% both fluorophores. Also plots the proportion of co-labeling as a function
% of the total number of MCs labeled.
%%
AREAs = unique([AniSliceTbl_avg.GFPSource;AniSliceTbl_avg.tdTomSource]);
SglSiteAvg = AniSliceTbl_avg(strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);

PropOverlap = cell(numel(AREAs),1);
for ii = 1:numel(AREAs)
    tDual = [SglSiteAvg(strcmp(SglSiteAvg.GFPSource,AREAs{ii})&strcmp(SglSiteAvg.tdTomSource,AREAs{ii}),:)];
    nLabeled{ii} = tDual.nGFP + tDual.nTdTom - tDual.nOverlap;
    PropOverlap{ii} = tDual.nOverlap ./ nLabeled{ii};    
end

figure; hold on
CLR = [70 130 180;255 125 0;255 200 0]./255;
FS = 20; MS = 16;
for ii = 1:numel(AREAs)
    subplot(121);hold on
    bar(ii,mean(PropOverlap{ii}),'facecolor',CLR(ii,:))
    plot(ii,PropOverlap{ii},'k.','markersize',MS)
    if ii == numel(AREAs)
        axis([0.25 3.75 0 1])
        set(gca,'xtick',1:3,'xticklabel',AREAs,'fontsize',FS*0.8)
        ylabel('Prop MCs co-labeled','fontsize',FS)
        xlabel('Single Site Mixed Virus Target','fontsize',FS)
    end
    subplot(122);hold on
    plot(nLabeled{ii},PropOverlap{ii},'.','color',CLR(ii,:),'markersize',MS*1.5)
    if ii == numel(AREAs)
        axis([0 325 0 1]);set(gca,'fontsize',MS)
        ylabel('Prop MCs co-labeled','fontsize',FS)
        xlabel('Number MCs Labeled','fontsize',FS)
        legend('BNST','MeA','PMCO','fontsize',FS)
    end
    
end