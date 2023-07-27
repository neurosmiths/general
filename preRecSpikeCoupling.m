%% [20170425] This script will examine spike phase coupling in epileptiform discharges trhoguhout the duration of the seizure. 
clear all
close all
set(0,'DefaultFigureRenderer','painters')

% direcvtory for saving figures. 
saveDir = '~/Figs/Seizures/spikeCoupling/';

% histograms
nBins = 29;
histEdges = -2*pi:(4*pi)./(nBins+1):2*pi;

%% first loading data - MUA and Broadband
nSeizures = 1; 
for sz = 1:nSeizures
	switch sz
		case 1
			patientID = 'C5-1';
			% load unit data
			unitDir = '/data/selected_data/sortedAPsMEA_emerix/c5_ictal_matched_units/'
			ptsz = 'c5_s1';
			unitList = dir([unitDir,ptsz,'*.mat'])
			recruited = 20;
		case 2
			patientID = 'C5-2';
		case 3
			patientID = 'C5-3';
		case 4
			patientID = 'NY442-1';
	end
	
	% using average LFP, ratehr than same channel. 	
	display('putting LFP in a matrix...')
	LFPmat = getLFPfromEdSrtSpks([unitDir,ptsz]);	
	[LFPbar] = DownSampleLFP(nanmean(LFPmat),3e4,2e3);

	% getting the LFP phase value for each unit time. 
	dF = dominantFrequency(patientID,LFPbar,2e3,[1 250],false);
	
	% looping over units
	nUnits = length(unitList);
	for un = 1:nUnits
		nUnits
		display(sprintf('calculating spike coupling for unit %d of %d ...',un,nUnits))
		% loading data fro each unit. 
		load([unitDir unitList(un).name])    	

		unitTimes = [unit.original.times; unit.matched.times];
		nSpikes = length(unitTimes);
		for sp = 1:nSpikes;
			% getting the value of each 
			[~,closestSample] = min(abs(repmat(unitTimes(sp),1,length(dF.tSec))-dF.tSec));
			
			% test condition and three controls. 
			freqCondition = 'test';
			switch freqCondition
				case {'test'}
					[~,freqBin] = min(abs(dF.fHz-repmat(dF.dominantFreq(closestSample),1,length(dF.fHz))));
				case {'gammaControl'}
					dF.fHz = 40;
					[~,freqBin] = min(abs(dF.fHz-repmat(dF.dominantFreq(closestSample),1,length(dF.fHz))));
				case {'HFAcontrol'}
					dF.fHz = 90;
					[~,freqBin] = min(abs(dF.fHz-repmat(dF.dominantFreq(closestSample),1,length(dF.fHz))));
			end
				
			if closestSample~=0
				tmp = squeeze(dF.PHIft);
				phaseVals{un}(sp) = tmp(freqBin,closestSample);
				if dF.tSec(closestSample)<recruited
					phaseVals_preRec{un}(sp) = tmp(freqBin,closestSample);
				elseif dF.tSec(closestSample)>=recruited
					phaseVals_postRec{un}(sp) = tmp(freqBin,closestSample);
				end
			end
		end
		
		%% for per-unit hotograms::
		% calculate histograms
		[phisto,edges] = histcounts(phaseVals{un},histEdges);
%		if isempty(phaseVals_postRec{un})
%			display(sprintf('no post-recruitment spikes for unit %d.',un))
%		else
			[phisto_preRec,edges_pR] = histcounts(phaseVals_preRec{un},histEdges);
			try
				[phisto_postRec,edges_Pr] = histcounts(phaseVals_postRec{un},histEdges);
			catch
				phisto_postRec = ones(1,100);
			end
			% [phisto_preTerm,edges] = histcounts(phaseVals_preTerm{un},histEdges);
	
			% plot histograms for each epoch.
			figure(un)
			rose(phaseVals_postRec{un},nBins+1);
			hold on
			rose(phaseVals_preRec{un},nBins+1);
			hold off
		
			% Testing for uniformity of distribution with Rayleigh's test. 
			pPre(un) = circ_rtest(phaseVals_preRec{un});
			pPost(un) = circ_rtest(phaseVals_postRec{un});

			% testing That the means are the same. [may be violated because distribution is bimodal] 
			pSameMean(un) = circ_wwtest(phaseVals_preRec{un},phaseVals_postRec{un});
			% testing that the distributions are the same. 
			pSameDist(un) = circ_ktest(phaseVals_preRec{un},phaseVals_postRec{un});

%			axis square tight
			maximize(un)
			fName = sprintf('%s_spikeCoupling.pdf',unitList(un).name(1:end-4));
			saveas(un,[saveDir fName])
			close(un)

 	%	end
		
	end
	
		%% for overall histogram:: TODO:: split this into pre and post recruitment epochs. 
	    display('now calculating unit coupling across the array.')
		[phisto_all_preRec,edges_all] = histcounts(cell2mat(phaseVals_preRec),histEdges);
    	[phisto_all_postRec,edges_all] = histcounts(cell2mat(phaseVals_postRec),histEdges);
    
        % plot histograms for each epoch.
        figure(9999)
		rose(cell2mat(phaseVals_postRec),nBins+1)
		hold on
  		rose(cell2mat(phaseVals_preRec),nBins+1)
        hold off

		% Testing for uniformity of distribution with Rayleigh's test. 
		pPre_all = circ_rtest(cell2mat(phaseVals_preRec))
		pPost_all = circ_rtest(cell2mat(phaseVals_postRec))
		
		% testing That the means are the same. [may be violated because distribution is bimodal] 
		pSameMean_all = circ_wwtest(cell2mat(phaseVals_preRec),cell2mat(phaseVals_postRec))
		% testing that the distributions are the same. 
		pSameDist_all = circ_ktest(cell2mat(phaseVals_preRec),cell2mat(phaseVals_postRec))

		% deets. 
        maximize(9999)
        fName = 'allUnits_spikeCoupling.pdf';
        saveas(9999,[saveDir fName])
        close(9999)

end

close all
figure(1)
hold on
histogram(pPre,100)       
histogram(pPost,100)       
hold off
saveas(1,[saveDir 'histogramOfPvalues.pdf'])




% So, took a while but new files are being uploaded now, in a directory called "c5_ictal_matched_units" in "/mnt/mfs/selected_data/sortedAPsMEA_emerix/". Each unit is a separate file which is pretty stupid but the result of the processing method. Might merge later. Anyway, each file has:

% details: just the usual details
% lfp: the unfiltered trace from that channel from time zero to +60 seconds, relative to all spike times in these files.
% unit: a structure of structures:
%     original: the original unit from proper spike sorting, with waveforms, times, and some metrics for reference.
%     matched: the waves and times from the ictal spikes that matched the original, plus metrics
%     possible: the waves and times that met the template matching criteria, but failed later tests, plus metrics.

% The metrics are line length, distance from centroid in principal component space, and correlation coefficient to the mean waveform from the original unit. In the "matched" and "possible" structures, it also includes line_length_probs and pc_space_probs, which are just the probability that that waveform corresponds to the original unit based on fitting a Gaussian distribution to the same metric in the original unit. This wasn't done for the correlation coefficient because it's non-Gaussian. All the raw metrics are in there too.

% Some of the template match results are a big mess, so I'm uploading another directory called figs within that one, showing the same figure as I shared previously for each unit on each channel. You can use that to easily discard the ones that clearly just matched garbage through the seizure, though I'm not 100% on what objective criteria you'd use for that, so maybe we can use the metrics to whittle down those messy channels to only the convincing waveforms. Nothing like a bit of subjective spike sorting...
