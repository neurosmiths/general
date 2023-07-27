function [MUAFile,microLFPFile,macroLFPfile] = preprocessSeizure(patientID, sz, fileName, filePath, Onset, Offset, preictalMins, postictalMins, inputFlag)
%PREPROCESSSEIZURE Processes seizures recorded on blackrock.
%   preprocessSeizures(patientID, sz, fileName, path, Onset, Offset,
%   preictalTime, PostictalTime) will output a preprocessed chunk of data
%   from seizure number [sz](scalar) from patient [patientID](string) with
%   [preictalTime](scalar) minutes of data before the [Onset](scalar) time
%   in samples and [postictalTime](scalar) minutes after the
%   [Offset](scalar) time in samples.
%
%   To save time, preprocessSeizure will ask you which preprocessed data
%   you would like to save. To avoid this step and process everything,
%   set inputFlag to 0;


% Author: Elliot H. Smith
% Version Date: 20160527
% https://github.com/elliothsmith/seizureCodes


% to Add:
%  - optional arguments for debugging MUA detection step
%  - spike sorting detected waveforms


%% some exceptions for a smoother UI
% fileName
if strcmp(fileName(end-4),'.')
    fileName = fileName(1:end-5);
end

% query user about preprocessing and saving options.
if inputFlag
    micromacroFlag = input('Preprocess 1.microelectrodes, 2.macro contacts, or 3.both?');
    if isequal(micromacroFlag,1)
        signalFlag = input('Preprocess 1.MUA, 2.LFP, or 3.both?');
        traceFlag = input('do you want to save downsampled MUA traces? (y/n)','s');
    elseif isequal(micromacroFlag,2)
        display('processing LFP on macro contacts');
        signalFlag = 3;
    elseif isequal(micromacroFlag,3)
        display('processing everything');
        signalFlag = 3;
        traceFlag = input('do you want to save downsampled MUA traces? (y/n)','s');        
    else
        error('Try again... \n\nPlease select 1, 2, or 3.')
    end
else
    micromacroFlag = 3;
    signalFlag = 3;
    traceFlag = 'y';
end


% depending on user response, preprocess macro and micro data.
if isequal(micromacroFlag,1)||isequal(micromacroFlag,3)
    %% This section preprocesses the microwire data.
    
    % loading data
    NS5 = openNSx([filePath fileName '.ns5']);
    Fs = 3e4;
    
    %% This section processses MUA:
    if isequal(signalFlag,1)||isequal(signalFlag,3)
        % detecting MUA
        MUA_BAND = [500 5000];
        [b,a] = fir1(250,MUA_BAND/(Fs/2));
        
        % loading data and checking to make sure that
        if size(NS5.Data,2)>Offset*15 + postictalMins*60*Fs
            display('not enough time at the end of this file to retain the requested postictal data.\nJust reading to end of file...')
            tmp = double(NS5.Data(:,Onset*15 - preictalMins*60*Fs:end));
            % timescale
            tsec = linspace(Onset*15 - preictalMins*60*Fs,length(tmp)./Fs,size(muamat,2));
        elseif (Onset*15 - preictalMins*60*Fs)<1
            tmp = double(NS5.Data(:,1:Offset*15 + postictalMins*60*Fs));
            display('not enough time at the beginning of this file to retain the requested postictal data.\nJust reading from the beginning of the file...')
        else
            tmp = double(NS5.Data(:,Onset*15 - preictalMins*60*Fs:Offset*15 + postictalMins*60*Fs));
        end
        
        % de-meaning and removing the common mode.
        dData = remove1stPC(tmp);
        clear tmp
        
        % filtering for MUA
        for ch = 1:size(NS5.Data,1)
            muamat(ch,:) = filtfilt(b,a,dData(ch,:));
            display(['filtering channel ' num2str(ch) ' of ' num2str(size(NS5.Data,1))])
        end
        
        % number of channels
        numChans = size(muamat,1);
        
        
        %% detect action potentials in the seizure files.
        DEBUG = 0;
        DETECTION_THRESHOLD = 5;
        ARTIFACT_THRESHOLD = 12;
        
        % filtering MUA trace
        WAVEFORM_RANGE = int32(floor(Fs*(-0.6)/1000)):int32(floor(Fs*1.0/1000));
        
        samplewin = floor(0.0005*Fs);
        avg_window = int16(60*Fs); % average threshold over this window
        
        % initialize the data structure
        mua_data.filter = MUA_BAND;
        mua_data.artifact_rejection_threshold = ARTIFACT_THRESHOLD;
        mua_data.detection_threshold = DETECTION_THRESHOLD;
        mua_data.thresholds = [DETECTION_THRESHOLD ARTIFACT_THRESHOLD];
        mua_data.start_time = Onset*15;
        mua_data.duration = size(muamat,2)/Fs;
        mua_data.nchannels = numChans;
        mua_data.nspikes = zeros(numChans,1);
        mua_data.fs = Fs;
        
        for c=1:numChans
            % calculate threshold based on interictal sample
            clear sig;
            sig = dData(c,1:avg_window);
            threshold = -DETECTION_THRESHOLD*std(sig);
            thresholds(c) = threshold;
            clear peaks;
            peaks = find_inflections (dData(c,:), 'minima');
            peaks = peaks(dData(c,peaks) < threshold);
            
            % remove peaks that are greater than X*SD of waveform max amplitude
            % greater than this is typically artifact
            maxabs = abs(mua(peaks));
            artifact_threshold = ARTIFACT_THRESHOLD*std(maxabs);
            
            %clear over;
            %over = find(maxabs > artifact_threshold);
            peaks(find(maxabs > artifact_threshold)) = [];
            
            % now save the waveforms
            waveforms = [];
            for k = 1:length(peaks)
                range = peaks(k) + WAVEFORM_RANGE;
                if min(range) < 1 || max(range) > length(mua)
                    % future:  for these, dip back into the original file
                    % otherwise, waveform analysis code can just skip
                    waveforms(k,:) = zeros(1, length(range));
                else
                    waveforms(k,:) = dData(c,range);
                end
            end
            mua_data.waveforms{c} = waveforms;
            
            % saving timestamps
            mua_data.timestamps{c} = peaks./Fs;
            mua_data.nspikes(c) = length(peaks);
            
            % debug feature.
            if DEBUG && length(peaks) > 0
                clf reset;
                subplot(3,1,1); hold on;
                title (['Channel ' num2str(c)])'
                x1 = 1:length(mua); x1 = x1./Fs;
                plot(x1,mua);
                for s = 1:length(peaks)
                    plot(x1(peaks(s)),threshold,'.r');
                end
                subplot(3,1,[2 3]); hold on;
                x2 = 1:49; x2 = x2./Fs; x2 = x2*1000;
                for s = 1:length(peaks)
                    plot(x2,waveforms(s,:));
                end
                title (['Threshold = ' num2str(threshold)]);
                pause;
            end
            
            % electrode labels
            trodeLabels{c} = strtrim(NS5.ElectrodesInfo.Label{c});
            
        end
        mua_data.thresholds = thresholds;
        
        
        %% [20160527] saving mua traces.
        if strcmp(traceFlag,'y')
            mua_data.muaTraces = dData(:,1:5:end);
            mua_data.timeSeconds = tsec(1:5:end);
        end
        
        
        %% saving MUA event times
        MUAFile = [patientID '_MUAtimes-' num2str(sz) '.mat'];
        save(MUAFile,'mua_data','trodeLabels','-v7.3')
        
        
    elseif isequal(signalFlag,2)||isequal(signalFlag,3)
        %% This section will preprocess the microelectrode LFP. \
        
       
%         if ~exist('dData','var')
%             % denoising full BW data by removing the first PC.
%             tmp = double(NS5.Data(:,Onset*15 - preictalMins*60*Fs:Offset*15 + postictalMins*60*Fs));
%             dData = remove1stPC(tmp);
%             clear tmp
%         end
%         
%         % downsample the seizure file and save to mat file with timestamps
%         Fs = 30000;
%         Fi = 1000;
%         display(['downasmpling LFP for seizure ' num2str(sz) ' from ' num2str(Fs) ' to ' num2str(Fi)])
%         DSFilter = fir1(256,Fi/(Fs/2));
%         
%         % initializing matrix
%         seizure_downFiltered = zeros(size(dData));
%         
%         % loooping over channels
%         for chs = 1:min(size(dData))
%             display(['Low pass filtering Channel ' num2str(chs) ' of 96'])
%             % noncausal fitlering
%             seizure_downFiltered(chs,:) = filtfilt(DSFilter,1,double(dData(chs,:)));
%         end
%         
%         % downsampling
%         seizure_downsampled = seizure_downFiltered(:,1:(Fs/Fi):end);
%         Seizure1K = seizure_downsampled;
%         
%         
%         %% [20150805] saving data from microwires.
%         microLFPFile = [patientID '_1Kdownsampled_microelectrodeLFP_seizure-' num2str(sz) '.mat'];
%         save(microLFPFile, 'Seizure1K', 'trodeLabels', '-v7.3')
        
        
    end
    
end

if isequal(micromacroFlag,2)||isequal(micromacroFlag,3)
    %% This section preprocesses the macro contacts
    clear HGpower Seizure2K
    
    %Fo
    Fo = 2e3;
    
    % update user
    display(['preprocessing LFP for macro contacts during seizure ' num2str(sz)])
    
    % load data
    NS3 = openNSx([filePath fileName '.ns3']);
    
    % snipping out requested segment
    tmp = double(NS3.Data(:,Onset - preictalMins*60*Fo:2:Offset + postictalMins*60*Fo));
    
    % timescale
    tsec = linspace( - preictalMins*60,((Offset/Fs) + postictalMins*60),size(Seizure1K,2));

    % demeaning data and removing the common mode. 
    Seizure1K = remove1stPC(tmp);
    clear tmp tmp1K
    Fs = 1e3;
    
    % High Gamma filter
    HG_BAND = [70 250];
    [b,a] = fir1(150,HG_BAND/(Fs/2));
    
    % loooping over channels
    for ch = 1:size(NS3.Data,1)
        HG = filtfilt(b,a,Seizure1K(ch,:));
        HGpower(ch,:) = abs(hilbert(HG));
        display(['filtering channel ' num2str(ch)])
        trodeLabels{ch} = NS3.ElectrodesInfo(ch).Label;
    end
    
    
    %% save processed data and elecatrode labels.
    macroLFPFile = [patientID '_2K_ECoGandHighGamma_sz' num2str(sz) '.mat'];
    save(macroLFPFile,'Seizure1K','HGpower','trodeLabels','tsec','-v7.3');
    
    
end


end


