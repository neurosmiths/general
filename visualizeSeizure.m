function [] = visualizeSeizure(patientID,szNum,Chan,saveFlag)
% VISUALIZESEIZURE plots seizure data from microelectrode recordings.
% 	visualizeSeizure(patientID,szNum,contact) searches the path for seizure szNum for patient patientID
%	and plots the data from macro contact number Chan and all microContacts.
%	to plot each and every macro contact, enter 'all' for Chan.
%
%	If saveFlag = 1, figures will be saved.
%
%   works great with data ouput from preprocessSeizure.m
%

% Author: Elliot H. Smith
% Version Date: 20160527
% https://github.com/elliothsmith/seizureCodes


%% checking for the existence of macroelectrode data.
if ischar(Chan) && strcmp(Chan,'all') && exist([patientID '_2K_ECoGandHighGamma_sz' num2str(szNum) '.mat'],'file')
    load([patientID '_2K_ECoGandHighGamma_sz' num2str(szNum) '.mat'])
    tsec = linspace(0,size(Seizure1K,2)./1e3,size(Seizure1K,2));
    %% plotting the seizure LFP and high gamma for each channel.
    for Chan = 1:size(Seizure1K,1)
        % plot LFP from macro contact in [Chan]
        figure(1)
        subplot(2,1,1)
        hold on
        plot(tsec,Seizure1K(Chan,:),'k');
        hold off
        ylabel('macroElectrode LFP (uV)', 'fontsize', 18)
        xlabel('time (seconds)', 'fontsize', 18)
        set(gca, 'linewidth', 2, 'fontsize', 16)
        title([strtrim(trodeLabels{Chan}) ' ||| seizure: ' num2str(szNum) ], 'fontsize', 18)
        axis tight
        
        % plot high gamma envelope on macro contact in [Chan]
        figure(1)
        subplot(2,1,2)
        hold on
        plot(tsec,HGpower(Chan,:),'k')
        ylabel('HighGamma Power (uV)', 'fontsize', 18)
        xlabel('time (seconds)', 'fontsize', 18)
        set(gca, 'linewidth', 2, 'fontsize', 16)
        axis tight
        
        % saving
        if saveFlag
            maximize(1)
            saveas(1,[patientID '_sz' num2str(szNum) '_macroContact' strtrim(char(trodeLabels{Chan})) '.pdf'])
        end
        
        close(1)
    end
    
elseif ~ischar(Chan) && exist([patientID '_2K_ECoGandHighGamma_sz' num2str(szNum) '.mat'],'file')
    load([patientID '_2K_ECoGandHighGamma_sz' num2str(szNum) '.mat'])
    
    % plot LFP from macro contact in [Chan]
    figure(1)
    subplot(2,1,1)
    hold on
    plot(tsec,Seizure1K(Chan,:))
    hold off
    ylabel('macroElectrode LFP (uV)', 'fontsize', 18)
    xlabel('time (seconds)', 'fontsize', 18)
    set(gca, 'linewidth', 2, 'fontsize', 16)
    title([strtrim(trodeLabels{Chan}) ' ||| seizure: ' num2str(szNum) ], 'fontsize', 18)
    axis tight
    
    % plot high gamma envelope on macro contact in [Chan]
    figure(1)
    subplot(2,1,2)
    hold on
    plot(tsec,HGpower(Chan,:),'k')
    ylabel('HighGamma Power (uV)', 'fontsize', 18)
    xlabel('time (seconds)', 'fontsize', 18)
    set(gca, 'linewidth', 2, 'fontsize', 16)
    axis tight
    
    % saving
    if saveFlag
        maximize(1)
        saveas(1,[patientID '_sz' num2str(szNum) '_macroContact' strtrim(char(trodeLabels{Chan})) '.pdf'])
    end
    
    close(1)
end


%% looking for multiunit activity data.
if isequal(exist([patientID '_MUAtimes-' num2str(szNum) '.mat'],'file'),2)
    % loading data.
    load([patientID '_MUAtimes-' num2str(szNum) '.mat'])
    
    % asking to plot MUA traces
    muaFlag = 'y'; % input('found MUA traces! would you like to visualize those? (y/n) \n     otherwise we"ll visualize the threshold crossings as a raster plot.','s')
    if strcmp(muaFlag,'y')
        
        % this section looks at the microelectrode labels to determine how
        % many BF electrodes there are.
        numBFs = floor(length(trodeLabels)/8);   % input('how many BFs were there in this pt?');
        for nbf = 1:numBFs
            % vertical spacing for MUA plots. 
            vertSpacing = 500;
            
            % plot broadband LFP from macro contact nearest to each micro
            % bundle.
            figure(nbf*10)
            % just starting with plotting the seizure.
            subplot(5,1,1)
            hold on
            plot(tsec,Seizure1K(((nbf-1)*8)+1,:)) % plots the medial electrode
            hold off
            ylabel('macroElectrode LFP (uV)', 'fontsize', 18)
            xlabel('time (seconds)', 'fontsize', 18)
            set(gca, 'linewidth', 2, 'fontsize', 16)
            title([strtrim(trodeLabels(((nbf-1)*8)+1)) ' ||| seizure: ' num2str(szNum) ], 'fontsize', 18)
            axis tight
            
            % plotting the micros.
            subplot(5,1,[2 3 4 5])
            hold on
            t30000 = linspace(0,size(mua_data.muaTraces,2)./3e4,size(mua_data.muaTraces,2));
            for cz = ((nbf-1)*8)+1:nbf*8
                plot(t30000,mua_data.muaTraces(cz,:)+(vertSpacing*cz),'k')
                text(0,vertSpacing*cz,trodeLabels(cz), 'fontsize', 16,'color',[0 1 0])
            end
            hold off
            ylabel('MUA traces')
            set(gca, 'linewidth', 2, 'fontsize', 16)
            axis tight
            ylim([((((nbf-1)*8)+1)*vertSpacing)-vertSpacing (vertSpacing*cz)+vertSpacing])
            
            % saving as a pdf and a .fig
            if saveFlag
                maximize(nbf*10)
                saveas(nbf*10,[patientID '_sz' num2str(szNum) '_microelectrodes' trodeLabels(((nbf-1)*8)+1) '.pdf'])
                saveas(nbf*10,[patientID '_sz' num2str(szNum) '_microelectrodes' trodeLabels(((nbf-1)*8)+1) '.fig'])
                close (nbf*10)
            end
            
            
        end
        
        % if not interested in MUA traces, plot the raster.
    elseif strcmp(muaFlag,'n')
        binRate = 1e3;
        
        % organize spike data.
        APs = zeros(length(trodeLabels),ceil(mua_data.duration.*binRate));
        for ch = 1:length(trodeLabels)
            APs(ch,ceil([mua_data.timestamps{ch}].*binRate)) = 1;
        end
        APs = logical(APs);
        
        tsec = linspace(1,ceil(mua_data.duration),size(APs,2));
        
        % plot raster
        figure(2)
        subplot(2,1,2)
        plotSpikeRaster(APs,'PlotType','vertline');
        hold on
        for cz = 1:length(trodeLabels)
            text(-10,cz,trodeLabels{cz}, 'fontsize', 16)
        end
        hold off
        set(gca, 'linewidth', 2, 'fontsize', 16)
        axis tight off
        
        % saving
        if saveFlag
            maximize(2)
            saveas(2,[patientID '_sz' num2str(szNum) '_microelectrodes.pdf'])
        end
        close(2)
        
    end
    
    %% TODO: add a figrue that plots a couple of micros from each bundle.
    
    
    if exist([patientID '_1Kdownsampled_microelectrodeLFP_seizure-' num2str(szNum) '.mat'],'file')
        load([patientID '_1Kdownsampled_microelectrodeLFP_seizure-' num2str(szNum) '.mat'])
        tsec = linspace(1,size(Seizure1K,2)./1e3,size(Seizure1K,2));
        
        % plot averaged LFP from each set of microelectrodes
        figure(2)
        subplot(2,1,1)
        plot(mean(Seizure1K),'k')
        xlabel('time (seconds)', 'fontsize', 18)
        ylabel('uElectrode LFP (uV)', 'fontsize', 18)
        set(gca, 'linewidth', 2, 'fontsize', 16)
        axis tight
        title([' seizure: ' num2str(szNum) ], 'fontsize', 18)
        
        % saving
        if saveFlag
            maximize(2)
            saveas(2,[patientID '_sz' num2str(szNum) '_microelectrodes.pdf'])
        end
        close(2)
    end
    
end


