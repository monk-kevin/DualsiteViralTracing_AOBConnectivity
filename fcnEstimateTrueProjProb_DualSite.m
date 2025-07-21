function fcnEstimateTrueProjProb_DualSite(AniSliceTbl_avg, nBoots)
%% Description
% This function takes a collection of projection data to estimate the true
% connectivity likelihood among pairs of areas when pairs of areas are
% infected with distinct fluorophores.
% 
% Specifically, for each % experimental sample, the observed co-labeled 
% rates are compared with ground-truth, probabilistic models to estimate 
% the range of true co-label % rates that is consistent with the observed 
% rates. For example, if a sample is composed of 100 labeled MCs, of which 
% 40 are co-labeled, then we can create probabilistic models of known 
% co-projection probabilities and determine for a sample of n = 100, which 
% ground-truth models are consistent with observing 40 co-labeled neurons. 
% The creation of these models are repeated nBoots (1000 by default) to 
% estimate the true underlying projection probabilities.
%
% These results are then plotted as a function of the number of neurons 
% labeled and color-coded by animal sex.
%% Input
% AniSliceTbl_avg: triplicate-averaged data for each experimental unit
% (animal/hemisphere pair). Created using the fcnCreateAniSliceTbl function
% nBoots: the number of unique times a probabilisitic model is created,
% 1000 by default.
%% Output
% A figure that shows for each pair of areas the range of ground-truth
% probabilities the observation is consistent with. Alternatively phrased,
% the vertical bars indicate the estimated underlying co-projection
% probabilities for neurons that project to multiple areas. Each row of
% panels shows the range for neurons that project to a given target area.
%%
figure

%definitions
AREAs = unique([AniSliceTbl_avg.GFPSource;AniSliceTbl_avg.tdTomSource]);
DualSiteAvg = AniSliceTbl_avg(~strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);
CLR = [70 130 180;255 125 0;255 200 0]./255;
FS = 20; MS = 16;
GroundTruthProp = 0:0.01:1;
GroundTruthPropHC = GroundTruthProp(1)-.005:.01:GroundTruthProp(end);

if nargin == 1
    nBoots = 1000;
end

nRow = 1;
LIMS = cell(numel(AREAs),numel(AREAs));
nLblCell = cell(numel(AREAs),numel(AREAs));
for ii = 1:numel(AREAs) % iterate across areas
    C = 0;
    for jj = 1:numel(AREAs) %iterate across partner areas
        if jj ~= ii
            C = C+1;
            tDual = [DualSiteAvg(strcmp(DualSiteAvg.GFPSource,AREAs{ii})&strcmp(DualSiteAvg.tdTomSource,AREAs{jj}),:);...
                DualSiteAvg(strcmp(DualSiteAvg.tdTomSource,AREAs{ii})&strcmp(DualSiteAvg.GFPSource,AREAs{jj}),:)];
            nCol = ceil(size(tDual,1)/nRow); 
            LIMS{ii,jj} = nan(size(tDual,1),2);
            AniSex{ii,jj} = tDual.AniSex;
            LimsVis = zeros(size(tDual,1),numel(GroundTruthProp));
            nLblCrossAni = nan(size(tDual,1),1);
            nLblAll = nan(size(tDual,1),1);
            for kk = 1:size(tDual,1)
                if strcmp(tDual.GFPSource{kk},AREAs{ii})
                    nLblAll(kk) = round(tDual.nGFP(kk));
                else
                    nLblAll(kk) = round(tDual.nTdTom(kk));
                end
            end
            [~,I] = sort(nLblAll);
            for KK = 1:size(tDual,1)
                kk = I(KK);
                nLbl = nLblAll(kk);
                propOverlap = tDual.nOverlap(kk) / nLbl;
                
                ChancePropLbl = nan(nBoots,numel(GroundTruthProp));
                for pp = 1:numel(GroundTruthProp)

                    for nn = 1:nBoots
                        RandLbl = rand(nLbl,1);
                        ChancePropLbl(nn,pp) = sum(RandLbl<=GroundTruthProp(pp))/nLbl;
                    end
                end
                PRTIL = [prctile(ChancePropLbl,5);prctile(ChancePropLbl,95)];
                LIMS{ii,jj}(kk,:) = [find(propOverlap>PRTIL(1,:)&propOverlap<PRTIL(2,:),1) find(propOverlap>PRTIL(1,:)&propOverlap<PRTIL(2,:),1,'last')];
                LimsVis(kk,LIMS{ii,jj}(kk,1):LIMS{ii,jj}(kk,2))=1;
                nLblCrossAni(kk) = nLbl;
            end
            [~,I] = sort(LIMS{ii,jj}(:,2));
            LIMS{ii,jj} = LIMS{ii,jj}(I,:);
            AniSex{ii,jj} = AniSex{ii,jj}(I);
            nLblCell{ii,jj} = nLblCrossAni(I);
            LimsVis = LimsVis(I,:);
            LimsVis_M = LimsVis(strcmp(AniSex{ii,jj},'M'),:);
            LimsVis_F = LimsVis(strcmp(AniSex{ii,jj},'F'),:);
            if C == 1
                Minds = 1;
                Finds = 7;
            else
                Minds = 4;
                Finds = 10;
            end
            
            % individual lines for each section triplicate, split by sex, 
            % ordered by number of labeled MCs and avg lines for upper/lower limits
            % log-scale
            
            subplot(6,6,[Minds+1 Finds+1] + 12*(ii-1))
            line(log(repmat(nLblCell{ii,jj}(strcmp(AniSex{ii,jj},'F')),1,2)'),...
                [LIMS{ii,jj}(strcmp(AniSex{ii,jj},'F'),1) LIMS{ii,jj}(strcmp(AniSex{ii,jj},'F'),2)]',...
                'color',CLR(jj,:)*0.75,'linewidth',4)
            line(log(repmat(nLblCell{ii,jj}(strcmp(AniSex{ii,jj},'M')),1,2)'),...
                [LIMS{ii,jj}(strcmp(AniSex{ii,jj},'M'),1) LIMS{ii,jj}(strcmp(AniSex{ii,jj},'M'),2)]',...
                'color',CLR(jj,:)*0.25,'linewidth',4)
            axis([0 6 0 100])
            ylabel('Ground Truth % Range','fontsize',FS)
            XTickVals = [1 5 10 25 50 100 200 300];
            set(gca,'fontsize',FS,'xtick',log(XTickVals),'XTickLabel',num2str(XTickVals'))
            xlabel('Number Labeled MCs','fontsize',FS)
            text(0.1,95,['AOB->' AREAs{ii} ' + AOB->' AREAs{jj}],'fontsize',FS,'color',CLR(jj,:),'fontweight','bold')
            
            subplot(6,6,[Minds+2 Finds+2] + 12*(ii-1))
            MnF = mean([LIMS{ii,jj}(strcmp(AniSex{ii,jj},'F'),1) LIMS{ii,jj}(strcmp(AniSex{ii,jj},'F'),2)]);
            MnM = mean([LIMS{ii,jj}(strcmp(AniSex{ii,jj},'M'),1) LIMS{ii,jj}(strcmp(AniSex{ii,jj},'M'),2)]);
            MnAll = mean([LIMS{ii,jj}(:,1) LIMS{ii,jj}(:,2)]);
            line([1 1],MnF,'color',CLR(jj,:)*0.75,'linewidth',4)
            line([1.5 1.5],MnM,'color',CLR(jj,:)*0.25,'linewidth',4)
            line([2 2],MnAll,'color',CLR(jj,:),'linewidth',8)
            text(0.55,115,['AvgRange_F: ' num2str(round(MnF(1))) '% - ' num2str(round(MnF(2))) '%'],'fontsize',FS,'color',CLR(jj,:)*0.75)
            text(0.55,105,['AvgRange_M: ' num2str(round(MnM(1))) '% - ' num2str(round(MnM(2))) '%'],'fontsize',FS,'color',CLR(jj,:)*0.25)
            text(0.55,95,['AvgRange_A_l_l: ' num2str(round(MnAll(1))) '% - ' num2str(round(MnAll(2))) '%'],'fontsize',FS,'color',CLR(jj,:))
            axis([0.5 2.5 0 120])
            ylabel('Avg Ground Truth % Range','fontsize',FS)
            set(gca,'xtick',1:.5:2,'XTickLabel',{'Female';'Male';'All'},'fontsize',FS,'ytick',0:20:100)

        end
    end
end
