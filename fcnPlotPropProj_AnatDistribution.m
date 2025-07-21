function fcnPlotPropProj_AnatDistribution(AniSliceTbl,AniSliceTbl_avg,SplitSexFlag)
%% Description
% Creates projection proportion figure that shows (1) the proportion of MCs
% that project to a given downstream area, (2) the anatomical distribution
% of labeled MCs as a function of the area they target [both in a
% normalized AOB map and in histogram] across the anterior-poster and
% dorsal-ventral axis

%% Inputs
% AniSliceTbl: non-averaged data table of labeled neurons
% AniSliceTbl_avg: averaged data table across triplicate sections
% both variables createed with the fcnCreateAniSliceTbl function
% SplitSexFlag: Binary flag that determines whether figure created will
% produce projections and anatomical distribution collapsed across sexes
% (SplitSexFlag = 0) or if figure produced will show proportion of
% projection MCs when the data is split across animal sex (SplitSexFlag =
% 1).
%% Outputs
% creates a figure with bar graphs, unit AOB plots, and histograms as
% described above.
% SplitSexFlag = 0: Proportions and anatomical distribution collapsed 
% across sex
% SplitSexFlag = 1: Proportion of projection MCs when splitting across male
% and female mice
%%
% divide table into dual and single site cohorts
DualSiteAvg = AniSliceTbl_avg(~strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);
SingleSiteAvg = AniSliceTbl_avg(strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);

% set up for iterations
AREAs = unique([DualSiteAvg.GFPSource;DualSiteAvg.tdTomSource]);
SEX = unique(DualSiteAvg.AniSex);
PropProj_sgle = cell(numel(AREAs),numel(SEX)+1);
PropProj_dual = cell(numel(AREAs),numel(SEX)+1);
PropProj_clps = cell(numel(AREAs),numel(SEX)+1);
PropProjMat_sgle = [];
PropProjMat_dual = [];
PropProjMat_clps = [];
nAni_sgle = nan(numel(AREAs),numel(SEX)+1);nHemi_sgle = nan(numel(AREAs),numel(SEX)+1);
nAni_dual = nan(numel(AREAs),numel(SEX)+1);nHemi_dual = nan(numel(AREAs),numel(SEX)+1);
nAni_clps = nan(numel(AREAs),numel(SEX)+1);nHemi_clps = nan(numel(AREAs),numel(SEX)+1);
AntPostHistProj = cell(numel(AREAs),1);

%iterate across individual areas to define the proportion of MCs that
%project downstream and the anatomical layout
for ii = 1:numel(AREAs)
    tDual = DualSiteAvg(contains(DualSiteAvg.GFPSource,AREAs{ii})|contains(DualSiteAvg.tdTomSource,AREAs{ii}),:);
    tSingle = SingleSiteAvg(contains(SingleSiteAvg.GFPSource,AREAs{ii})&contains(SingleSiteAvg.tdTomSource,AREAs{ii}),:);
    %proportion of projections split by sex
    for jj = 1:numel(SEX)+1
        % split by sex
        if jj < numel(SEX)+1
            t1Dual = tDual(strcmp(tDual.AniSex,SEX{jj}),:);
            t1Single = tSingle(strcmp(tSingle.AniSex,SEX{jj}),:);
        else
            %collapsed across animal sex
            t1Dual = tDual;
            t1Single = tSingle;
        end

        % define the total number labeled neurons as nTd + nGFP - nOverlap
        propTotalLabeled = t1Single.propTdTom+t1Single.propGFP-t1Single.propOverlap;

        PropProj_sgle{ii,jj} = propTotalLabeled;
        PropProj_dual{ii,jj} = [t1Dual.propGFP(strcmp(t1Dual.GFPSource,AREAs{ii})); t1Dual.propTdTom(strcmp(t1Dual.tdTomSource,AREAs{ii}))];
        PropProj_clps{ii,jj} = [t1Dual.propGFP(strcmp(t1Dual.GFPSource,AREAs{ii})); t1Dual.propTdTom(strcmp(t1Dual.tdTomSource,AREAs{ii}));...
            propTotalLabeled];

        tempMat_sgle = PropProj_sgle{ii,jj};
        tempMat_sgle(:,2) = ii;
        tempMat_sgle(:,3) = jj;
        PropProjMat_sgle = [PropProjMat_sgle;tempMat_sgle];
        nAni_sgle(ii,jj) = numel(unique(t1Single.AniID));
        nHemi_sgle(ii,jj) = size(t1Single,1);

        tempMat_dual = PropProj_dual{ii,jj};
        tempMat_dual(:,2) = ii;
        tempMat_dual(:,3) = jj;
        PropProjMat_dual = [PropProjMat_dual;tempMat_dual];
        nAni_dual(ii,jj) = numel(unique(t1Dual.AniID));
        nHemi_dual(ii,jj) = size(t1Dual,1);

        tempMat_clps = PropProj_clps{ii,jj};
        tempMat_clps(:,2) = ii;
        tempMat_clps(:,3) = jj;
        PropProjMat_clps = [PropProjMat_clps;tempMat_clps];
        nAni_clps(ii,jj) = nAni_sgle(ii,jj)+nAni_dual(ii,jj);
        nHemi_clps(ii,jj) = nHemi_sgle(ii,jj)+nHemi_dual(ii,jj);
    end

    %combine anatomical data across sections for population analyses
    AntPostHistProj{ii} = [tDual.AntPostHC_GFP(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.AntPostHC_tdTom(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
        tSingle.AntPostHC_GFP;tSingle.AntPostHC_tdTom];

    AntPostHist_NT{ii} = [tDual.AntPostHC_NT(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.AntPostHC_NT(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
        tSingle.AntPostHC_NT;tSingle.AntPostHC_NT];

    AntPostHistProj_norm{ii} = AntPostHistProj{ii} - AntPostHist_NT{ii};

    DorsVentHistProj{ii} = [tDual.DorsVentHC_GFP(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.DorsVentHC_tdTom(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
        tSingle.DorsVentHC_GFP;tSingle.DorsVentHC_tdTom];

    DorsVentHist_NT{ii} = [tDual.DorsVentHC_NT(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.DorsVentHC_NT(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
        tSingle.DorsVentHC_NT;tSingle.DorsVentHC_NT];

    DorsVentHistProj_norm{ii} = DorsVentHistProj{ii} - DorsVentHist_NT{ii};

end

figure
MS = 24;FS = 16;
CLR = [70 130 180;255 125 0;255 200 0]./255;
if SplitSexFlag == 0
    % plot the proportion of MCs that project to a given area
    jj = 3; % plotting data for collapsed across sex
    subplot(141)
    hold on
    StatsMat = [];
    for ii = 1:numel(AREAs)
        bar(ii,cellfun(@mean,PropProj_clps(ii,jj)),'facecolor',CLR(ii,:))
        plot(ones(size(PropProj_clps{ii,jj}))*ii,PropProj_clps{ii,jj},'k.','markersize',MS)
        tmp = [PropProj_clps{ii,jj} ones(size(PropProj_clps{ii,jj}))*ii];
        StatsMat = [StatsMat;tmp];
    end
    errorbar(cellfun(@mean,PropProj_clps(:,3)),cellfun(@std,PropProj_clps(:,3))./sqrt(cellfun(@numel,PropProj_clps(:,3))),'color','k','linewidth',6,'linestyle','none')
    axis([0.5 3.5 0 .75])
    set(gca,'xtick',1:3,'xticklabel',AREAs,'fontsize',FS)
    [p,~,stats] = anova1(StatsMat(:,1),StatsMat(:,2),'off');

    % posthoc tests comparing distributions of projection probabilities
    c = multcompare(stats,'Display','off');
    line([1 2],[.65 .65],'color','k','linewidth',2)
    line([1 3],[.7 .7],'color','k','linewidth',2)
    text(1,0.675,['*: p = ' num2str(round(c(1,end),3))],'fontsize',FS)
    text(1,0.725,['**: p = ' num2str(round(c(2,end),3))],'fontsize',FS)
    ylabel('Proportion of MCs','fontsize',FS)
    xlabel('Target Area','fontsize',FS)

    % heatmaps of AOB locations on Unit AOB
    img_dim = [81 27];
    AntPostDorsVentMatGFP = nan(img_dim(1),img_dim(2),size(AniSliceTbl,1));
    AntPostDorsVentMatTdTom = nan(img_dim(1),img_dim(2),size(AniSliceTbl,1));
    for ii = 1:size(AniSliceTbl,1)
        x = AniSliceTbl.FracAntPostGFP{ii}(:,1);y = AniSliceTbl.FracAntPostGFP{ii}(:,2);
        img = zeros(img_dim);

        xt = floor(rescale(x, 1, img_dim(2)));
        yt = floor(rescale(y, 1, img_dim(2)));
        idx = sub2ind(img_dim,yt,xt); %assumes square img_dim (i.e., img_dim(1) = img_dim(2))
        img(idx) = 1;
        AntPostDorsVentMatGFP(:,:,ii) = flipud(img);

        x = AniSliceTbl.FracAntPostTdTom{ii}(:,1);y = AniSliceTbl.FracAntPostTdTom{ii}(:,2);
        img = zeros(img_dim);

        xt = floor(rescale(x, 1, img_dim(2)));
        yt = floor(rescale(y, 1, img_dim(2)));
        idx = sub2ind(img_dim,yt,xt);
        img(idx) = 1;
        AntPostDorsVentMatTdTom(:,:,ii) = flipud(img);

    end

    % plotting unit AOB to show anatomical heatmaps of AOB MCs
    CLIM = [0 .66];
    for ii = 1:numel(AREAs)
        MNmat = cat(3,...
            mean(AntPostDorsVentMatGFP(:,:,strcmp(AniSliceTbl.GFPSource,AREAs{ii})),3),...
            mean(AntPostDorsVentMatTdTom(:,:,strcmp(AniSliceTbl.GFPSource,AREAs{ii})),3));
        subplot(7,4,ii+1)
        imagesc(1:img_dim(1),1:img_dim(2),mean(MNmat,3),CLIM)
        title(['Area: ' AREAs{ii}])
        set(gca,'xtick',[1 img_dim(1)],'ytick',[1 img_dim(2)],'xticklabel',{'P','A'},'yticklabel',{'D','V'},'fontsize',FS)
    end

    % histogram of AP axis
    % subplot(2,6,10:12); hold on
    nSteps = size(AntPostHistProj{1},2);

    FaceAlpha = 1;
    for ii = 1:numel(AREAs)
        % subplot(4,6,ii+12);hold on
        subplot(7,4,[5 9 13] + ii);hold on
        plot(mean(AntPostHist_NT{ii}),'color','k','linewidth',2)
        plot(mean(AntPostHistProj{ii}),'color',CLR(ii,:),'linewidth',2)


        fill([1:nSteps nSteps:-1:1],...
            [mean(AntPostHist_NT{ii}) + (std(AntPostHist_NT{ii})./sqrt(size(AntPostHist_NT{ii},1)))...
            fliplr(mean(AntPostHist_NT{ii}) - (std(AntPostHist_NT{ii})./sqrt(size(AntPostHist_NT{ii},1))))],...
            'k','edgecolor','none','facealpha',FaceAlpha)
        fill([1:nSteps nSteps:-1:1],...
            [mean(AntPostHistProj{ii}) + (std(AntPostHistProj{ii})./sqrt(size(AntPostHistProj{ii},1))) ...
            fliplr(mean(AntPostHistProj{ii}) - (std(AntPostHistProj{ii})./sqrt(size(AntPostHistProj{ii},1))))],...
            CLR(ii,:),'edgecolor','none','facealpha',FaceAlpha)

        axis([0.5 nSteps+.5 0 .2])
        set(gca,'xtick',[1 nSteps],'xticklabel',{'P';'A'},'fontsize',FS-2)
        title(['Area: ' AREAs{ii}],'fontsize',FS)

        % subplot(4,6,ii + 18);hold on
        subplot(7,4,[17 21 25] + ii);hold on
        plot(mean(DorsVentHist_NT{ii}),1:nSteps,'color','k','linewidth',2)
        plot(mean(DorsVentHistProj{ii}),1:nSteps,'color',CLR(ii,:),'linewidth',2)
        fill([mean(DorsVentHist_NT{ii}) + (std(DorsVentHist_NT{ii})./sqrt(size(DorsVentHist_NT{ii},1))) ...
            fliplr(mean(DorsVentHist_NT{ii}) - (std(DorsVentHist_NT{ii})./sqrt(size(DorsVentHist_NT{ii},1))))],...
            [1:nSteps nSteps:-1:1],...
            'k','edgecolor','none','facealpha',FaceAlpha)
        fill([mean(DorsVentHistProj{ii}) + (std(DorsVentHistProj{ii})./sqrt(size(DorsVentHistProj{ii},1))) ...
            fliplr(mean(DorsVentHistProj{ii}) - (std(DorsVentHistProj{ii})./sqrt(size(DorsVentHistProj{ii},1))))],...
            [1:nSteps nSteps:-1:1],CLR(ii,:),...
            'edgecolor','none','facealpha',FaceAlpha)
        axis([0 .2 0.5 nSteps+.5])
        set(gca,'ytick',[1 nSteps],'yticklabel',{'V';'D'},'fontsize',FS-2)
    end

    % stats for anatomical distribution (comparing labeled vs neurotrace)
    for ii = 1:numel(AREAs)
        for jj = 1:nSteps
            pRS_AP(ii,jj) = ranksum(AntPostHistProj{ii}(:,jj),AntPostHist_NT{ii}(:,jj));
            pRS_DV(ii,jj) = ranksum(DorsVentHistProj{ii}(:,jj),DorsVentHist_NT{ii}(:,jj));

            pDiffFrom0_AP(ii,jj) = signrank(AntPostHistProj_norm{ii}(:,jj));
            pDiffFrom0_DV(ii,jj) = signrank(DorsVentHistProj_norm{ii}(:,jj));
        end
    end

    ALPHA = .05/nSteps;
    Yval = [0.175 0.175 0.175];
    Xval = [0.175 0.175 0.175];

    for ii = 1:numel(AREAs)
        % rank sum
        Xinds = find(pRS_AP(ii,:)<ALPHA);
        Yinds = find(pRS_DV(ii,:)<ALPHA);

        subplot(7,4,[5 9 13] + ii);hold on
        if ~isempty(Xinds)
            plot(Xinds,Yval(ii),'.','markersize',MS,'color',CLR(ii,:))
        end
        legend('All MCs','Proj. MCs')

        subplot(7,4,[17 21 25] + ii);hold on
        if ~isempty(Yinds)
            plot(Xval(ii),Yinds,'.','markersize',MS,'color',CLR(ii,:))
        end

    end

else % split across sex
    ClrScale = [1 0.65];
    AntPostHistProj_sex = cell(numel(AREAs),1);
    for ii = 1:numel(AREAs)
        tDual = DualSiteAvg(contains(DualSiteAvg.GFPSource,AREAs{ii})|contains(DualSiteAvg.tdTomSource,AREAs{ii}),:);
        tSingle = SingleSiteAvg(contains(SingleSiteAvg.GFPSource,AREAs{ii})&contains(SingleSiteAvg.tdTomSource,AREAs{ii}),:);

        AntPostHistProj{ii} = [tDual.AntPostHC_GFP(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.AntPostHC_tdTom(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
            tSingle.AntPostHC_GFP;tSingle.AntPostHC_tdTom];

        AntPostHistProj_sex{ii} = [tDual.AniSex(strcmp(tDual.GFPSource,AREAs{ii}),:);tDual.AniSex(strcmp(tDual.tdTomSource,AREAs{ii}),:);...
            tSingle.AniSex;tSingle.AniSex];

        for jj = 1:2 %for the two sexes
            subplot(1,3,ii)
            hold on
            bar(jj,mean(PropProj_clps{ii,jj}),'facecolor',CLR(ii,:)*ClrScale(jj))
            plot(ones(size(PropProj_clps{ii,jj}))*jj,PropProj_clps{ii,jj},'k.','markersize',MS)
        end

        subplot(1,3,ii)
        errorbar(cellfun(@mean,PropProj_clps(ii,1:2)),cellfun(@std,PropProj_clps(ii,1:2))./sqrt(cellfun(@numel,PropProj_clps(ii,1:2))),'color','k','linewidth',2,'linestyle','none')
        axis([0.5 2.5 0 0.85])
        [~,p_t2]=ttest2(PropProj_clps{ii,1},PropProj_clps{ii,2});
        line([1 2],[0.7 0.7],'color','k','linewidth',2)
        text(1.15,0.75,['p = ' num2str(round(p_t2,3))],'fontsize',FS)
        ylabel({AREAs{ii};'Proportion MCs'})
        set(gca,'xtick',1:2,'xticklabel',{'Female','Male'},'fontsize',FS)
    end
end