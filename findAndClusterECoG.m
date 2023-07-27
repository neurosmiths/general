% let's leave this as a script for now.
clear all
close all

set(0,'DefaultFigureRenderer','painters')


%% first loading data - MUA and Broadband
nSeizures = 4;
for sz = 2:nSeizures
    % clearing data:: save everything in the ptData struct to be safe!
    clearvars -except ptData nSeizures sz
    
    switch sz
        % Utah Array Pts.
            
    end
    
    
    %% detecting discharges or loading previously-detected discharges.
    dDfile = sprintf('~/Data/Seizures/detectedDischarges/%s_detectedDischarges.mat',patientID);
    if ~exist(dDfile,'file')
        [dD] = detectDischarges(patientID,data,Fs,[],coreChan); % don't rename the output variable or you have to do this again!!
    else
        load(dDfile)
        dD = dischargeData;
    end
    
    
    %% finding minimum time window around discharges.
    IDIwin = round(min(diff(dD.times)),3);
    IDIsamps = ceil((IDIwin./2)*Fs);
    nDischarges = length(dD.times);
    
    
    %% sorting channels by the amount of high gamma.
    nChans = length(chanLabels);
    
    
    %% need to do all spectral analyses here particulalry need to rank HFA channels.
    % my plan is to find the dominnant frequencies and then rank the channels by the amount of HFA following recruitment. Using the dominant frequencies might be a useful for determining recruitment times. [TODO subsequently test how these features generalize to new ECoG data sets.- This could be a second paper -]
    [dF] = dominantFrequency(patientID,data,Fs,[],false);
    dFname = sprintf('~/Data/Seizures/%s_dominantFrequencyData.mat',patientID);
    save(dFname,'-v7.3')
    
    
    %% determining HFA values at discharge times and ranking based on postRecValues.
    tD.dFampsAtDischargeTimes = dF.HFA(:,~isnan(dD.timeIndices));
    tD.preRecDischargAmps = dF.HFA(:,~isnan(dD.timeIndices(dD.tSec<recruited)));
    tD.coreAmps = dF.HFA(:,~isnan(dD.timeIndices(dD.tSec>recruited)));
    
    [~,tD.channelRanks] = sort(nanmean(tD.coreAmps,2),'ascend');
    % TODO:: rank these for both pre and post recruitment epochs
    % [should be a cool pattern!!!
    % [also required for tSNE
    
    
    % making a dope colormap.
    cMap = hsv2rgb([repmat(linspace(0,5/6,nDischarges)',nChans,1) repmat(linspace(1,0,nChans)',nDischarges,2)]);
    % TODO: make Saturation/Value a function of distance from ictal core.
    
    
    %% building a voltage tensor, X: [channels x data x discharges] and a spectral tensor, S: [channels X frequency X samples X discharges]
    X = [];
    S = [];
    for ds = 1:nDischarges
        [~,closestSample] = min(abs(repmat(dD.times(ds),1,length(dD.tSec))-dD.tSec));
        samps = closestSample-IDIsamps:closestSample+IDIsamps;
        X(:,:,ds) = data(tD.channelRanks,samps); % X:[rankedChannels X samples X discharges]
        S(:,:,:,ds) = dF.Sft(tD.channelRanks,:,samps); % S: [channels X frequency X samples X discharges]
        tX(ds,:) = dD.tSec(samps); % tX [time X discharges]
    end
    Xmat = tmatn(X,2);
    vecLabels = repmat(chanLabels(tD.channelRanks),1,nDischarges)';
    labelNums = repmat(tD.channelRanks<=nChans,nDischarges,1)';
    
    
    %% entropy
    entropy = true;
    if entropy
        nBins = 50;
        for d2 = 1:nDischarges
            updateUser(d2,50,nDischarges,'calculating entropy for discharge')
            for c2 = 1:nChans
                updateUser(c2,10,nChans,'calculating entropy for channel')
                % differential
                Hd(d2,c2) = differentialEntropy(X(c2,:,d2),nBins);
                
                % spectral
                tmpS = squeeze(nanmean(S(c2,dF.fHz<150,:,d2),3))./(1./dF.fHz(dF.fHz<150));
                [Hs(d2,c2),Sbins] = spectralEntropy(tmpS);
                clear tmpS
            end
        end
        
        % plotting
        figure(22)
        colormap redblue
        plotmultipleaxes(1,1,2,0.03,22)
        imagesc(Hd')
        axis tight xy
        colorbar
        set(gca,'yticklabel',chanLabels,'ytick',1:nChans)
        
        plotmultipleaxes(2,1,2,0.03,22)
        imagesc(Hs')
        axis tight xy 
        colorbar
        set(gca,'yticklabel',chanLabels,'ytick',1:nChans)
        
        ylabel('electrodes')
        xlabel('discharges')
        
        % saving entropy matrices for quantification across patients. 
        ptData(sz).Hd = Hd;
        ptData(sz).Hs = Hs;
        
        % saving figures.
        maximize(22)
        pause(2)
        saveas(22,sprintf('~/Figs/Seizures/entropy/%s_spectralanddifferentialEntropy.pdf',patientID))
        pause(2)
        close(22)
    end
    
        
    
    %     % ripe for cool visualization of discharges...
    %     figure(1)
    %     subplot(2,1,1)
    %     colormap(cMap)
    %     plot(linspace(-IDIwin/Fs,IDIwin/Fs,length(samps)),Xmat)
    %     set(gca,'Xticklabel',{})
    %     ylabel('discharge amplitude (uV)')
    %     axis tight
    
    
    tSNE = false;
    if tSNE
        %% clustering discharges (w/out feature selection) using tSNE.
        tD.init_dims = 50;
        tD.perplexity = 30;
        tSName = sprintf('~/Data/Seizures/tSNE_discharges/%s_2DtSNEmapping_allECoGdischarges__perp%d__initDims%d.mat',patientID,tD.perplexity,tD.init_dims);
        if ~exist(tSName,'file')
            tD.dateTimeStarted = datetime;
            tD.cMap = cMap;
            tD.labelNums = [];
            tD.labels = vecLabels;
            tD.no_dims = 2;
            
            tic
            if ~isempty(tD.labelNums)
                figure
                halfMaximize(gcf)
            end
            tD.mappedXmat = tsne(Xmat',tD.labelNums,tD.no_dims,tD.init_dims,tD.perplexity);
            tD.dateTimeCompleted = datetime;
            save(tSName)
            A = toc;
            display(sprintf('t-SNE clustering took %s seconds.',A))
        else
            load(tSName)
        end
        
        
        %% plotting clusters
        cMap = hsv2rgb([repmat(linspace(0,5/6,nDischarges)',nChans,1) repmat(linspace(1,0,nChans)',nDischarges,2)]);
        figure(1)
        %     subplot(2,1,2)
        scatter(tD.mappedXmat(:,1),tD.mappedXmat(:,2),10,cMap,'filled')
        axis off tight
        halfMaximize(1,'left')
        
        saveas(1,sprintf('~/Figs/Seizures/tSNE/%s_tSNE_allDischargesandChannels__perp%d__initDims%d.pdf',patientID,tD.perplexity,tD.init_dims))
    end
    
    
    %% plot discharges from each cluster including spectra
    % TODO:: make a dope color plot of stacked discharges or HFA
    pltDischarges = false;
    savePath = '~/Figs/Seizures/detectedDischarges/';
    perChan = true; %true plots all channels. False plots means.
    if pltDischarges
        for dd = 1:nDischarges
            if perChan
                for cc = 1:nChans
                    border = 0.05;
                    figgyH = cc;
                    figure(figgyH)
                    ah1 = plotmultipleaxes(1,1,2,border,figgyH);
                    plot(tX(dd,:),X(cc,:,dd),'k')             % X:[rankedChannels X samples X discharges]
                    xlabel('time (s)')
                    ylabel('LFP (uV)')
                    axis tight
                    
                    % plotting spectra.
                    ah2 = plotmultipleaxes(2,2,2,border,figgyH);
                    imagesc(tX(dd,:),dF.fHz(dF.fHz<150),squeeze(S(cc,dF.fHz<150,:,dd))./repmat(1./dF.fHz(dF.fHz<150)',1,length(tX(dd,:))))  % S: [channels X frequency X samples X discharges]
                    colormap jet
                    colorbar
                    xlabel('time (s)')
                    ylabel('frequency (Hz)')
                    zlabel('normalized power')
                    axis square xy
                    
                    ah3 = plotmultipleaxes(4,2,2,border,figgyH);
                    hold on
                    %                     patch([squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3)) fliplr(squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3)))],[dF.fHz(dF.fHz<150) fliplr(dF.fHz(dF.fHz<150))],rgb('salmon'),'facealpha',0.5,'edgecolor','none')
                    patch([(squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3))./(1./dF.fHz(dF.fHz<150)))+(squeeze(nanstd(S(cc,dF.fHz<150,:,dd),[],3))./(1./dF.fHz(dF.fHz<150)))...
                        fliplr((squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3)./(1./dF.fHz(dF.fHz<150))))-(squeeze(nanstd(S(cc,dF.fHz<150,:,dd),[],3))./(1./dF.fHz(dF.fHz<150))))]...
                        ,[dF.fHz(dF.fHz<150) fliplr(dF.fHz(dF.fHz<150))],rgb('salmon'),'facealpha',0.5,'edgecolor','none')
                    %                     plot(squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3)),dF.fHz(dF.fHz<150),'color',rgb('salmon'),'linewidth',2)
                    plot(squeeze(nanmean(S(cc,dF.fHz<150,:,dd),3))./(1./dF.fHz(dF.fHz<150)),dF.fHz(dF.fHz<150),'color',rgb('salmon'),'linewidth',2)
                    hold off
                    ylabel('frequency (Hz)')
                    xlabel('normalized power')
                    axis square
                    
                    %saving
                    maximize(figgyH)
                    fName = sprintf('%s_channel%d_discharge%d_VoltageSpectra.pdf',patientID,cc,dd);
                    saveas(figgyH,[savePath fName]);
                    close(figgyH)
                end
            end
        end
    end
    %%
    
    
end


