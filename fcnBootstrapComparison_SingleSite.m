function fcnBootstrapComparison_SingleSite(AniSliceTbl_avg,nBoots,AlphaPrctile)
%% Description
% This function compares the observed co-labeling proportions in single-
% site animals with what would be expected if neurons were independently
% targeting downstream areas. Note that this is done on a sample-by-sample
% basis. 
% 
% Namely, for each sample for each area, we create a dummy population of
% neurons with a size defined by the total number of units observed to
% project to that area for a given sample. We then find the number of dummy
% neurons that are co-labeled by random chance. This process is repeated a
% number of times (nBoots, 1000 by default). This distribution tells us the
% proportion of neurons expected to be co-labeled and we can compare the
% whether the observed proportion of co-labeling exceeds this expectation.
%
% For visualization purposes, bootstrapped distributions are normalized
% with Z-score normalization, but statistical tests are performed on
% non-normalized data.
%% Inputs
% AniSliceTbl_avg: triplicate-averaged data for each experimental unit
% (animal/hemisphere pair). Created using the fcnCreateAniSliceTbl function
%% Output
% Figure demonstrating the proportion of times that an observed
% distribution exceeds co-projection probabilities as defined through
% bootstrap distributions.
%%
AREAs = unique([AniSliceTbl_avg.GFPSource;AniSliceTbl_avg.tdTomSource]);
SingleSiteAvg = AniSliceTbl_avg(strcmp(AniSliceTbl_avg.GFPSource,AniSliceTbl_avg.tdTomSource),:);
if nargin == 1
    nBoots = 1000;
elseif nargin == 2
    AlphaPrctile = 99;
end

figure
nPairs = numel(AREAs)*2;
PropObsExceedingThr = nan(nPairs,1);
PRS = [1 2;1 3;2 3];
CC = 0;
HistLim = [-2 10.5];HistBin = 0.5;
CLR = [70 130 180;255 125 0;255 200 0]./255;
FS = 16; MS = 24;
ObsSampleAll = cell(0);
for ii = 1:numel(AREAs)
    C = 0;
    ObsSampleAll{ii} = [];
    for jj = 1:2
            C = C+1;
            CC = CC+1;
            tSingle = SingleSiteAvg(strcmp(SingleSiteAvg.GFPSource,AREAs{ii}),:);
            
            ChancePropLbl_Z = nan(nBoots,size(tSingle,1));
            Obs_Z = nan(size(tSingle,1),2);
            ObsSample = nan(size(tSingle,1),2);
    
            for kk = 1:size(tSingle,1)
                ChancePropLbl = nan(nBoots,1);
                
                % create bootstrap distriubtion
                for nn = 1:nBoots
                    if jj == 1
                        RandLbl = rand(round(tSingle.nGFP(kk)),1);
                        ChancePropLbl(nn) = sum(RandLbl<=tSingle.propTdTom(kk))/numel(RandLbl); % the proportion of green neurons expected to also be red based on the proportion of red neurons
                    else
                        RandLbl = rand(round(tSingle.nTdTom(kk)),1);
                        ChancePropLbl(nn) = sum(RandLbl<=tSingle.propGFP(kk))/numel(RandLbl); % the proportion of red neurons expected to also be green based on the proportion of green neurons
                    end
                end
                ChancePropLbl_Z(:,kk) = (ChancePropLbl - mean(ChancePropLbl)) ./ std(ChancePropLbl); %z-normalizing chance distribution

                % normalizing observed distributions
                if jj == 1
                    Obs_Z(kk,jj) = (tSingle.nOverlap(kk)/tSingle.nGFP(kk) - mean(ChancePropLbl)) ./ std(ChancePropLbl);
                    ObsSample(kk,:) = [tSingle.nOverlap(kk)/tSingle.nGFP(kk) tSingle.nGFP(kk)];
                else
                    Obs_Z(kk,jj) = (tSingle.nOverlap(kk)/tSingle.nTdTom(kk) - mean(ChancePropLbl)) ./ std(ChancePropLbl);
                    ObsSample(kk,:) = [tSingle.nOverlap(kk)/tSingle.nTdTom(kk) tSingle.nTdTom(kk)];
                end
            end


            ObsSampleAll{ii} = [ObsSampleAll{ii};ObsSample];
            subplot(3,3,C + 3*(ii-1))
            hold on
            histogram(mean(ChancePropLbl_Z,2),HistLim(1):HistBin:HistLim(2),'Normalization','probability','facecolor',ones(1,3)*.6,'facealpha',1)
            % plotting fit line over gray histograms
            muData = mean(mean(ChancePropLbl_Z,2));
            stdData = std(mean(ChancePropLbl_Z,2));
            binCenters = HistLim(1):HistBin/10:HistLim(2);
            y = normpdf(binCenters,muData,stdData);
            % pd = fitdist(mean(ChancePropLbl_Z,2),'normal');
            % y = pdf(pd,binCenters);
            PDscale = max(histcounts(mean(ChancePropLbl_Z,2),HistLim(1):HistBin:HistLim(2),'Normalization','probability'));
            % paramci(pd,'alpha',.01);
            % plot(binCenters,y*PDscale,'k','linewidth',2)
            % plot(binCenters,y,'k','linewidth',2)

            % xline(prctile(mean(ChancePropLbl_Z,2),99),'linestyle','--','color','r','linewidth',2)
            histogram(Obs_Z(:,jj),HistLim(1):HistBin:HistLim(2),'Normalization','probability','facecolor',CLR(jj,:),'facealpha',1)
            if jj == 1
                title({['AOB-->' AREAs{ii} ' that also'];['-->' AREAs{ii} ': GFP']},'fontsize',FS)
            else
                title({['AOB-->' AREAs{ii} ' that also'];['-->' AREAs{ii} ': tdTom']},'fontsize',FS)
            end

            axis([-3 11 0 1])
            xlabel('Co-Projection Probability (Z-Norm to Bootstrap)')

            PropObsExceedingThr(CC) = sum(Obs_Z(:,jj)' >= prctile(ChancePropLbl_Z,AlphaPrctile)) / numel(Obs_Z(:,jj));

            % z-scored
            coeff = polyfit(Obs_Z(:,jj),ObsSample(:,2),1);
            xFit = linspace(min(Obs_Z(:,jj)),max(Obs_Z(:,jj)));
            yFit = polyval(coeff,xFit);
            [r_corr,p_corr] = corrcoef(Obs_Z(:,jj),ObsSample(:,2));
            % subplot(3,4,3 + 4*(ii-1))
            % hold on
            % plot(Obs_Z(:,jj),ObsSample(:,2),'.','markersize',MS,'color',CLR(jj,:))
            % plot(xFit,yFit,'color',CLR(jj,:),'linewidth',2)
            % 
            % if C == 1
            %     text(0.1, 275, ['R^2 = ' num2str(round(r_corr(2)^2,2)) ', p = ' num2str(round(p_corr(2),2))],'color',CLR(jj,:)*0.67)
            % else
            %     text(0.1, 250, ['R^2 = ' num2str(round(r_corr(2)^2,2)) ', p = ' num2str(round(p_corr(2),2))],'color',CLR(jj,:)*0.67)
            % 
            %     axis([-2 11 0 300])
            %     xlabel('Observed Prop Overlap (Znorm)','fontsize',FS)
            %     ylabel('Number of Back-Labeled Neurons','fontsize',FS)
            % end
    end
end

subplot(1,3,3);hold on
ClrInd = [1 1 2 2 3 3];
XInd = [1 2 4 5 7 8];
C = 0;
for ii = 1:CC
    bar(XInd(ii),PropObsExceedingThr(ii),'FaceColor',CLR(ClrInd(ii),:))
    if rem(ii,2)~=0
        C = C+1;
        text(XInd(ii),1.05,AREAs{C},'fontsize',FS)
    end
end
set(gca,'xtick',XInd,'XTickLabel',{'gfp','tdtom','gfp','tdtom','gfp','tdtom'},'fontsize',FS,'XTickLabelRotation',45,'ylim',[0 1.1],'ytick',0:.1:1)
xline([3 6],'color','k','linestyle','--','linewidth',2)
xlabel('fluorophore')
ylabel(['Proportion of Sections that Exceed ' num2str(AlphaPrctile) 'th Percentile'])