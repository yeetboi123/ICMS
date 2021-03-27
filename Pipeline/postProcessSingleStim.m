%% First load the firings data

addpath basic_util/ mdaio/ neuroshare/ 
parent_dir = uigetdir('Date folder')
input_dir = fullfile(parent_dir,'Processing')
image_save_path = fullfile(parent_dir, 'output_images')
addpath(fullfile(parent_dir,'raster_helper_functions/'))

BIN_RESOLUTION = 0.005; % time resolution of histogram bin in seconds
PRESTIM_BLANK_TIME= 0.000 ; % zero out output t seconds before stim
POSTSTIM_BLANK_TIME = 0.000; % zero out output t seconds after stim
POSTSTIM_AVG_WINDOW = 0.005; % window of post-stim average waveform 
PRE_STIM_WINDOW = 0.15; % window before and after stim for plotting rasterplot 
POST_STIM_WINDOW = 0.15; % window before and after stim for plotting rasterplot 

fs = 3e4;

%% Load input files

addpath(input_dir); 

load('channelInfo.mat') % contains info about each recording channel
load('stimChannelInfo.mat'); % contains info about each STIM channel 
load('sorted_TS_current.mat') % contains info about the stim trials 
    
A = readmda('firings.curated_4_5.mda'); % read firings data
Ch = A(1,:); 
timeStamp_samples = A(2,:); timeStamp_seconds = timeStamp_samples./fs;
neuron = A(3,:);
clear A; 

% cluster_info includes information from recording channels 
cluster_info = zeros(1,3);
idx1= unique(neuron);
NEURON_NUM = numel(idx1);

STIM_CHANNEL_NUM = numel(unique(sorted_TS_current(:,3)));

for i = 1:numel(idx1)
    neuron_idx = idx1(i);
    cluster_info(i,1) = neuron_idx; % not important 
    primaryCh = Ch(find((neuron==neuron_idx),1));
    cluster_info(i,2) = primaryCh; % map this to original channel name 
    cluster_info(i,3) = numel(find(cluster_info(:,2)==primaryCh));
end

filtered_data=readmda(strcat(input_dir,'/filt.mda'));

stim_times = sorted_TS_current(:,1);

%% Create cell array that contains stim pulse timestamps for every stim channel/current combination and replicant 

getCurrentInd = @(x) find(abs(sorted_TS_current(:,2) - x) < 0.001);
getChInd = @(x) find(sorted_TS_current(:,3) == x);
getReplicantTS = @(stim_ch,stim_current) sorted_TS_current(intersect(getCurrentInd(stim_current),getChInd(stim_ch)),1);
new_order = unique(sorted_TS_current(:,3));

% Find the replicants for each stimulation channel 
% Find the timestamps for neuron i for each replicant 
% Order the replicants so that first index is recording channel, second is
% current, third is replicant. NOTE: INFO IS FROM STIMULATING CHANNELS ONLY!!!
nCH = STIM_CHANNEL_NUM; nCU = 4; currents = [0.5,1,2,5]; NUM_REPLICANTS = 10;
stim_pulse_arr = cell(nCH,nCU);
for i = 1:nCH
    for j = 1:nCU
            stim_ch = new_order(i); stim_current = currents(j);
            stim_pulse_arr{i,j} = getReplicantTS(stim_ch,stim_current);
    end
end

%% Process each cluster or putative neuron
h = figure(1); clf;
h.Position = [61 194 1278 1084];
stimulating_ch = {};
stim_graph = zeros(nCH,nCH);
% for neuron_idx = 1:NEURON_NUM
for neuron_idx = 1:NEURON_NUM
    
    primaryCh = cluster_info(neuron_idx,2);
    IDforCh = cluster_info(neuron_idx,3); % there can be more than 1 cluster/neuron for each primary channel 
    
    neuron_idx_info = channelInfo(primaryCh,:); 
    name = neuron_idx_info{1}; ch_ind = neuron_idx_info{2}; x_coord = neuron_idx_info{3}; y_coord = neuron_idx_info{4};
    
    figure(1); clf;
    t = sgtitle(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh)),'Interpreter','none');
    set(t, 'horizontalAlignment', 'right')
   
    % selectBoolean selects all indices in which the unit is equal to the unit of interest
    selectBoolean = (neuron == idx1(neuron_idx));
    % find the primary channel that this unit was recorded on using the
    % first instance in which selectBoolean equals true
    primaryChannel = Ch(find((selectBoolean==1),1));
    % unitTimestamps holds the timestamps of the unit of interest
    neuronTimestamps = timeStamp_seconds(selectBoolean)';
    neuronTimestampsSamples = timeStamp_samples (selectBoolean)';
    % each cell in group will hold the timestamps of the units within the
    % period before and after each stimulation (in other words, the first
    % cell in group will hold all unit timestamps for the window
    % before the stim pulse and the window after stimulation)
    % group is a 1xT cell for T trials/stimulus pulses

    % Key helper function!
    % Keeps firings if timestamps fall in user-defined time boundaries
    [FIRING_DATA ,post_stim_avg_window_1_array,post_stim_avg_window_2_array] = ...
            getEvents2(neuronTimestamps, stim_pulse_arr, POSTSTIM_BLANK_TIME, POSTSTIM_AVG_WINDOW, ...
            PRE_STIM_WINDOW,POST_STIM_WINDOW,nCH,nCU,NUM_REPLICANTS,currents, new_order);
        
    % plotting the raster plot on left-hand side
    % each row in the rasterplot is the spikes in a stimulation window
    generateRasterplot2(FIRING_DATA, stimChannelInfo, PRE_STIM_WINDOW,POST_STIM_WINDOW,nCH,nCU,NUM_REPLICANTS,0,[]);

    numPoints = (PRE_STIM_WINDOW + POST_STIM_WINDOW)/BIN_RESOLUTION ;
    time_hist = linspace(-PRE_STIM_WINDOW,POST_STIM_WINDOW,numPoints+1);
    subplot_tight(1,4,3)
    colors = get(gca,'colororder');

    for fig = 1:nCH
%         subplot_tight(nCH,4,(nCH*4-3)-4*(fig-1)+2,[0.01,0.01]);
        
        combinedPSTH = [];
        curr_firings = [];
        
        for curr = 1:nCU

            % take average of replicants

            for rep = 1:NUM_REPLICANTS
                temp = [FIRING_DATA{fig,curr}{rep,1}{1,2}];
                curr_firings = [curr_firings; temp];
            end    

            [allCounts,~,event_binInd] = histcounts(curr_firings,time_hist);
           [preStim_count, postStim_count, counts,preStimMean,preStimStd] = getHistLimits(numPoints, allCounts,time_hist);
            
            postStim_limit = floor(0.02/BIN_RESOLUTION);
            postStim_end = postStim_count+postStim_limit;
            if floor(postStim_end) ~= postStim_end
                postStim_end = ceil(postStim_end);
            end
            postStim = counts(postStim_count:postStim_end);

            pre_temp = counts(1:preStim_count);
            post_temp = counts(postStim_count:end);
            PSTH = [pre_temp post_temp];
            combinedPSTH  = [combinedPSTH PSTH];
            
            rec_ch = primaryCh; 
            stim_ch = fig;
            
            if preStimMean*post_temp(1) > 3*preStimStd || preStimMean*post_temp(2) > 3*preStimStd
                stim_graph(stim_ch,rec_ch) = stim_graph(stim_ch,rec_ch)+curr;
            else
                stim_graph(stim_ch,rec_ch) = 0;
            end

        end

        plot((30*combinedPSTH)+120 * (fig-1),'LineWidth',1,'color',colors(7,:))
        
        
        hold on;
        gridxy(preStim_count + ((length(time_hist)-1)* [0:3]),'LineStyle',':','linewidth',0.25);
        hold on;
        title('Peristimulus Time Histogram');
        xlabel('Time(s)');
        ylabel('Normalized firing rate');
        set(gca,'visible','off')
        axis tight
    end

    %% plot the average waveform of unit in a small inset at the top right
    % hard coded # of samples before and after firing timestamps
    tBefore = 75; tAfter = 75;
    generateAverageWaveform(neuronTimestampsSamples,...
        filtered_data,primaryCh,tBefore,tAfter);
% 
% 
%     %% plot the average waveform of "specific" bin in a small inset at the top right
%     % the histogram outputs the histogram bin index FOR ALL UNITS WITHIN
%     % THE LIMITS OF THE STIMULATION. Therefore, we need to pare down the
%     % detected unit timestamps to just those within 5 msec of the
%     % stimulation times.
% 
    generateAvgPostStimWaveform(post_stim_avg_window_1_array,post_stim_avg_window_2_array,...
                         filtered_data,primaryCh,tBefore,tAfter)

    % keep array of relevant units with 2nd column being physical channels
    % and 3rd column being the index for multiple units for 1 channel
%         relevant_unit = [unitNum, cluster_info(unitNum,3)];
%         relevant_units = [relevant_units; relevant_unit];
    % name each image with relative channel and index of unit for relative
    % channel
    cd(image_save_path);
%     saveas(h,char(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh), '.png')))
    export_fig(char(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh), '.png')),'-native')
    cd('..');

end
cd('..');