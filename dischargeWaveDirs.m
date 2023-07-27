function [ArrayOrder,propAng,propMag,propDirection] = dischargeWaveDirs(data,times,tsec,map,patient,burstPlot,MOVIE)

% DISCHARGEWAVEDIRS plots the direction vectors of epilepstiform discharges
%   across a grid of electrodes.
%
%   inputs: data = channels x samples data matrix
%           times = 1 x n times in seconds at which to look for a wave
%           tsec = a vector represetning time in seconds (length = samples)
%           map = a map (matrix) of electrode numbers (a la electrodepinout.m)
%
%   optional inputs:
%           'patient' = patient identifier. 
%           'burstPlot' = true or false (1/0; 0 default) for plotting the data
%           'Movie' = true/false (1/0; 0 defalut) for making a movie
%
%   outputs:ArrayOrder = ordered electrodes per burst
%           propAng = angle of velocity vectors (polar)
%           propMag = magnitude of velocity vectors (polar) 
%           propDirection = direction and magnitude of the vector in
%               cartesian coordinates. (compass plot) 


% Author: Elliot H. Smith
% Version Date: 20160518


if ~exist('burstPlot','var')
    burstPlot = 0;
end

if ~exist('MOVIE','var')
    MOVIE = 0;
end

if ~exist('patient','var')
    patient = num2str(date);
end

Fs = 3e4; % round(length(tsec>0)./tsec(end));

% how wide is the discharge (in samples)?
burstHalfWidth = 50./Fs;

% making a colormap for bursts
burstCmap = colormap(summer(length(times)));


%% making a movie of propagation directions across the grid to compare with microelectrode array
% setting up movie
if MOVIE
    movie_file_name = ['./' patient '_travelingWaveDirections.avi'];
    movie_object = VideoWriter(movie_file_name,'Motion JPEG AVI');
    movie_object.FrameRate = 5;
    movie_object.Quality = 50;
    open(movie_object);
    
    % setting up figure
    movie_handle = figure(999);
    figure(999)
    maximize(999)
    
    % plot seizure
    subplot(2,2,[1 2])
    plot(tsec, mean(data),'k','linewidth',1);
    axis off
    
    % set up compass plot
    subplot(2,2,4)
    h0 = compass(0,2);
    set(h0,'Color',rgb('white'));
    axis manual
    hold on
end


% looping over bursts.
for bst = 1:length(times)
    
    bindices = tsec>times(bst)-burstHalfWidth & tsec<times(bst)+burstHalfWidth ;
    [~,mindices] = min(diff(data(:,bindices)'));  % mindices: the indices of the minima for each channel
    
    [~,channelOrder] = sort(mindices,'ascend');    % channelOrder: channels ordered by their burst minimae
    channelBurstOrder(:,bst) = channelOrder;
    
    
    %% %%%%%%%%%%%%%%%% Visualizing Burst Slopes %%%%%%%%%%%%%%%%
    Cmap = colormap(cool(size(data,1)));
    
    %% converting burst orders into a map.
    Ich = zeros(length(channelBurstOrder(:,bst)));
    Jch = zeros(length(channelBurstOrder(:,bst)));
    
    % looping over the bursts
    for chz = 1:size(data,1)
        try
        [I,J] = ind2sub(size(map),find(map==channelBurstOrder(chz,bst)));
        % saving subscripts over Channels
        Ich(chz) = I;
        Jch(chz) = J;
        
        ArrayOrder(I,J,bst) = chz;
        catch 
           display(sprintf('couldn"t find channel %d in the map.',chz)); 
        end
        
        % plotting waves across ECoG electrodes.
        if (burstPlot && bst>length(times)-30)
            figure(bst)
            subplot(2,1,1)
            hold on
            plot(data(channelBurstOrder(chz,bst),bindices)','color',Cmap(chz,:));
            scatter(5+(J*2),-200+(I*40),200,'MarkerEdgeColor','none','MarkerFaceColor',Cmap(chz,:));
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% calculating burst direction.
    binaryImage = ArrayOrder(:,:,bst)~=0;
    measurementsPlus = regionprops(logical(binaryImage),ArrayOrder(:,:,bst),'WeightedCentroid');
    measurementsMinus = regionprops(logical(binaryImage),size(data,1)-ArrayOrder(:,:,bst),'WeightedCentroid');
    
    Rminus = measurementsMinus.WeightedCentroid;
    Rplus = measurementsPlus.WeightedCentroid;
    
    propDirection(:,bst) = Rplus-Rminus;
    
    if MOVIE;
        % indicating which burst
        figure(999)
        subplot(2,2,[1 2])
        hold on
        plot(times(bst),0,'color',burstCmap(bst,:),'marker','*')
        hold off
        
        % plotting delay map
        subplot(2,2,3)
        imagesc(mean(ArrayOrder(:,:,bst),3))
        
        % plotting the burst directions
        subplot(2,2,4)
        hE = compass(propDirection(1,bst),propDirection(2,bst));
        hE.Color = burstCmap(bst,:);
        hE.LineWidth = 2;
        
        % writing video
        writeVideo(movie_object,getframe(movie_handle));
    end
    
end


% angles in polar coordinates
[propAng,propMag] = cart2pol(propDirection(1,:),propDirection(2,:));

