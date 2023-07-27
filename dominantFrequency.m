function [dominantFreq] = dominantFrequency(patientID,data,Fs,fPass,plotFlag)
% DOMINANTFREQUENCY determines frequency peaks of a singal from wavelet
% decomposition.
%
%	[dominantFreq] = dominantFrequency(patientID,data,Fs,fPass,plotFlag)
%
% patientID is a string identifying a patient.
% data is an array of data (i.e. [channels X samples])
% Fs is the sampling frequency.
% fPass is a 2 element vector specifying the frequency window to be examined. default: [1 Fs/2]
% plotFlag is a boolean flag indicating whether to plot

% Author:: [EHS20170602]

dominantFreq.tSec = linspace(0,length(data)/Fs,length(data));

% default window. 
if isempty(fPass)
	fPass = [1 Fs/2];
end

% calculating wavelet spectrogram.
nChans = size(data,1);
for ch = 1:nChans
    display(sprintf('spectral calculations for channel %d out of %d',ch,nChans))
    [wcoeffs,~,scale] = basewave4(data(ch,:),Fs,fPass(1),fPass(2),8,0);
    
    % determining scales in Hz
    scaleFreqs = linspace(fPass(1),fPass(2),length(scale));
    %scaleScales = scale*Fs/2;
    dominantFreq.fHz = scaleFreqs;
    
    % power and related transformations here. Sft:[chan X freq X time]
    dominantFreq.Sft(ch,:,:) = squeeze(abs(hilbert(wcoeffs)));
	dominantFreq.PHIft(ch,:,:) = squeeze(angle(hilbert(wcoeffs)));
    SftTilt = squeeze(dominantFreq.Sft(ch,:,:)./repmat(1./scaleFreqs,[1,1,length(dominantFreq.Sft(ch,:,:))]));
    
    % do dominant frequency stuff here.
    %     [~,maxY] = max(squeeze(Sft(ch,scaleFreqs<fPass(2)-50,:)),[],2);
    [~,maxX] = max(squeeze(dominantFreq.Sft(ch,scaleFreqs<fPass(2)-50,:)));
    [~,maxXtilt] = max(squeeze(SftTilt(scaleFreqs<fPass(2)-50,:)));
    fMat = repmat(scaleFreqs,1,length(dominantFreq.Sft));
    dominantFreq.dominantFreqRaw(ch,:) = fMat(maxX);
    dominantFreq.dominantFreqTilt(ch,:) = fMat(maxXtilt);
    dominantFreq.SmoF = Fs/4;
    dominantFreq.spectralTilt = '1/f';
    dominantFreq.dominantFreq(ch,:) = smooth(fMat(maxX),dominantFreq.SmoF);
    
    % define HFA band and rank channels?
    tmpHFAwin = [50 200];
    dominantFreq.HFAwin = tmpHFAwin;
    dominantFreq.HFA(ch,:) = squeeze(nanmean(dominantFreq.Sft(ch,scaleFreqs>tmpHFAwin(1) & scaleFreqs<=tmpHFAwin(2),:),2)); % HFA:[chans x time];
    
    if plotFlag
        border = 0.05;
        % plotting the raw voltage and high gamma
        figure(1123)
        ah1 = plotmultipleaxes(1,1,2,border,1123);
        hold on
        plot(dominantFreq.tSec,data(ch,:),'color',rgb('gray'),'linewidth',1)
        plot(dominantFreq.tSec,dominantFreq.HFA(ch,:),'color',rgb('darkgreen'),'linewidth',1)
        line([0 0],[-1000 0],'linewidth',5,'color',rgb('rosybrown'))
        text(0,-500,'1 mV')
        ylabel('LFP (uV)','fontsize',18)
        set(gca,'linewidth',2,'fontsize',16)
        title([patientID ' channel: ' num2str(ch)])
        axis off
        hold off
        
        % plotting the wavelet decomposition
        ah2 = plotmultipleaxes(2,1,2,border,1123);
        hold on
        imagesc(dominantFreq.tSec, scaleFreqs(scaleFreqs<fPass(2)-50), squeeze(dominantFreq.Sft(ch, scaleFreqs<fPass(2)-50, :)),[0 max(max(squeeze(dominantFreq.Sft(ch,:,:))))])
        plot(dominantFreq.tSec,dominantFreq.dominantFreq(ch,:))
        colormap hot
        colorbar
        zlabel('power (dB)','fontsize',16)
        ylabel('frequency (Hz)','fontsize',16)
        xlabel('time (seconds)','fontsize',16)
        ylim([1 fPass(2)-50])
        set(gca,'linewidth',2,'fontsize',16)
        axis xy tight
        hold off
        
        % saving
        savePath = '~/Figs/Seizures/';
        savestr = sprintf('%s_channel_%d_ECoGspectrogram.pdf',patientID,ch);
        
        maximize(1123)
        
        saveas(1123,[savePath savestr])
        close(1123)
        
    end
end

dominantFrequencyPlot = false;
if dominantFrequencyPlot
    figure(1124)
    surf(dominantFreq.tSec,1:nChans,dominantFreq.dominantFreq)
    
    ylabel('frequency (Hz)','fontsize',16)
    xlabel('time (seconds)','fontsize',16)
    ylim([1 fPass(2)-50])
    set(gca,'linewidth',2,'fontsize',16)
    axis xy tight
    
    title(sprintf('%s dominant frequencies',patientID));
    
    % saving
    savePath = '~/Figs/Seizures/';
    savestr = sprintf('%s_dominantFreqs_allChans.pdf',patientID,ch);
    
    
    maximize(1124)
    
    saveas(1124,[savePath savestr])
    close(1124)
end
