function [dischargeData] = detectDischarges(patientID,data,Fs,method,channel,thresholdScaling)
% DETECTDISCHARGES detects discharges in seizures using LFPs.
%
%	[dischargeData] = detectDischarges(patientID,data,Fs,method,channel);
%
%   patientID is a string identifying the patient.
% 	data is a matrix of multichannel electrophysiology data [channels X samples(@Fs samples per second)]
%	method is a string indicating which method to use to detect discharges.
%
%       now:
%           'cluster'
%           'moving'
%
%   signal determines which data to operate on.
% 
%	now:
% 			'high gamma'
% 			'MUA'
% 			'firing rate'
% 			'local minima'
% 			'inflections'
% 			'predetermined times'
%
% 	channel specifies the channel to use to find discharge times.
%   
%   thresholdScaling is a scalar between 0 and 1 that can adjust the
%   threshold proportionately. 
% 
%	This function outputs a structure with discharge times, amplitudes, overall time, and the time/date that the data were processed. 
% 
%  

% probably have TODO some more stuff here.

if ~exist('thresholdScaling','var')
    thresholdScaling = 1;
end

% define filter specs
% fWin = [3 20];% works pretty well using low frequencies.
fWin = [50 200];
b = fir1(150,fWin./(Fs/2));

% filtering all of the channels simultaneously.
nChans = size(data,1);
for ch = 1:nChans
    % updating user
    display(sprintf('filtering channel %d of %d between %d and %d Hz.',ch,nChans,fWin(1),fWin(2)));
    HFA(ch,:) = abs(hilbert(filtfilt(b,1,data(ch,:))));
end

% detecting from a single channel or averaged data?
if isempty(channel)
    HFAbar = mean(HFA);
else
    HFAbar = HFA(channel,:);
end

% dischrge detection parameters.
if fWin(1)<40
    eta = 0.5; 		  % amplitude trheshold in percentage of std of signal.
elseif fWin(1)>=40
    eta = 2;
end
smoothType = 'Gaussian';
if strcmp(smoothType,'moving')
    smool = 50;
    HFAbarSmo = smooth(HFAbar,'moving',smool);
elseif strcmp(smoothType,'Gaussian')
    sigmar = 10; % gaussian width in samples. 
    Gauss = fspecial('gaussian',[1 8*sigmar],sigmar);
    HFAbarSmo = conv(HFAbar,Gauss,'same');
end

% detecting dischrges.
tSec = linspace(0,size(data,2)./Fs,size(data,2));

% first using amplitude...
[~,HGpeaks] = findpeaks(HFAbarSmo);
peakIdcs = false(1,length(HFAbarSmo));
peakIdcs(HGpeaks') = true;
HGenvelope = NaN(1,length(HFAbarSmo));
HGenvelope(peakIdcs) = HFAbarSmo(peakIdcs);
ampThresh = (thresholdScaling)*eta*nanstd(HGenvelope);
HGenvelope(HGenvelope<ampThresh) = NaN;

% saving discharge times and amplitudes.
dischargeData.tSec = tSec;
dischargeData.timeIndices = HGenvelope;
dischargeData.times = tSec(~isnan(HGenvelope));
dischargeData.amplitudes = HGenvelope(~isnan(HGenvelope));

% then pruning discharges based on timing.
% first removing outliers.
longDurations = outliers(diff(dischargeData.times));
figure(1)
subplot(2,1,1)
halfMaximize(1,'right')
hold on
plot(tSec,HFAbar,'k')
plot(tSec,HFAbarSmo,'color',rgb('darkseagreen'))
line([tSec(1) tSec(end)],[ampThresh ampThresh])
plot(tSec,HGenvelope,'*r')
plot(dischargeData.times(longDurations),dischargeData.amplitudes(longDurations),'*b')
hold off
xlabel('time (seconds)')
ylabel('signal (black: raw // sea green: smoothed)')
title('remove temporally outlying discharges (blue)?')
firstDischarges = input('remove the first N discharges?')
lastDischarges = input('remove the last N discharges?')

% removing early and late outliers.
dischargeData.times(1:firstDischarges) = [];
dischargeData.amplitudes(1:firstDischarges) = [];
dischargeData.times(end-lastDischarges:end) = [];
dischargeData.amplitudes(end-lastDischarges:end) = [];

% [20170530] tried a gaussian mixture model to dissect the bimodality of
% the distribution, but just finding local minima in a smoothed histogram
% works just as well, empirically.
figure(1)
subplot(2,1,2)
hold on
scatter([0 diff(dischargeData.times)],[0 diff(dischargeData.amplitudes)],10,'filled')
X = 0:0.01:0.4;
[N,Xhat] = histcounts(diff(dischargeData.times),X);
plot(X,[0 N*10])
xlabel('inter-discharge interval (s)')
ylabel('interdischarge amplitude difference (uV)')

% clustering based temporal segregation.
% TODO:: implement sliding window approach.
discardVector = [];
while isempty(discardVector)
    K = 7; % number of clusters (i recommend an odd number. you'll see why)
    XX = [[0 diff(dischargeData.times)]; [0 diff(dischargeData.amplitudes)]];
    [Ls,~] = kmeans(XX,K);
    cMap = lines(K);
    subplot(2,1,1)
    plot(tSec,HFAbar,'k')
    for p = 1:length(Ls)
        subplot(2,1,1)
        hold on
        scatter(dischargeData.times(p),dischargeData.amplitudes(p),20,cMap(Ls(p),:),'filled')
        hold off
        subplot(2,1,2)
        hold on
        scatter(XX(1,p),XX(2,p),12,cMap(Ls(p),:),'filled')
        hold off
    end
    for lg = 1:K
        text(0.4,(K-lg)*100,num2str(lg),'color',cMap(lg,:))
    end
    text(-100,0.3,'cluster numbers color-coded above this text','color',[0 0 0])
    title('Detections colored by kluster. Which to discard? [] => re-cluster')
    halfMaximize(1,'right')
    % ask user which clusters to discard
    discardVector = input('which clusters should we discard ([] to recluster, 99 to accept result)?');
    
    if ~isequal(discardVector,99)
        for dc = 1:length(discardVector)
            dischargeData.times(Ls==discardVector(dc)) = [];
            dischargeData.amplitudes(Ls==discardVector(dc)) = [];
            Ls(Ls==discardVector(dc)) = [];
        end
        
        % PLOT FINAL RESULT
        close(1)
        figure(1)
        hold on
        plot(tSec,HFAbar,'k')
        scatter(dischargeData.times,dischargeData.amplitudes,30,rgb('deeppink'),'filled')
        hold off
        xlabel('time (s)')
        ylabel('signal amplitude')
        halfMaximize(1,'right')
        if ~isempty(discardVector)
            title('do you want to re-cluster? [Y/n]')
            reCluster = input('do you want to re-cluster? [Y/n]','s');
            if strcmp(reCluster,'Y')
                discardVector=[];
            elseif strcmp(reCluster,'n')
                display('accepting result...')
            else
                display('please enter "Y" or "n"')
                discardVector=[];
            end
        end
    end
end

dischargeData.interDischargeInterval = [0 diff(dischargeData.times)];
dischargeData.acceptedTime = datetime;

% % sliding window based
% [~,I] = max(N);
% winWidth = X(I+2);

saveFlag = input('result accepted! /nsave figure in ~/Figs/Seizures/ and data in ~/Data/Seizures? [Y/n]','s');
if strcmp(saveFlag,'Y')
    saveas(gcf,sprintf('~/Figs/Seizures/%s_dischargeDetectionResult.pdf',patientID))
    fName = sprintf('%s_detectedDischarges.mat',patientID);
    save(fullfile('/mnt/mfs/selected_data/elliotWorking/Data/Seizures/detectedDischarges_new/',fName),'dischargeData')
end


