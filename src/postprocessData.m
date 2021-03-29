function []=postprocessData()

    addpath utils/
    addpath mdaio/ neuroshare/ raster_helper_functions/

    processing_path = uigetdir('Select directory to read and write files for postprocessing');
    image_save_path = fullfile(processing_path,'output_images');

    BIN_RESOLUTION = 0.005; % time resolution of histogram bin in seconds
    PRESTIM_BLANK_TIME= 0.005 ; % zero out output t seconds before stim
    POSTSTIM_BLANK_TIME = 0.000; % zero out output t seconds after stim
    POSTSTIM_AVG_WINDOW = 0.005; % window of post-stim average waveform 
    PRE_STIM_WINDOW = 0.1; % window before and after stim for plotting rasterplot 
    POST_STIM_WINDOW = 0.1; % window before and after stim for plotting rasterplot 

    addpath(processing_path)
    
    channelInfo=load('channelInfo.mat').channelInfo;
    stimChannelInfo = channelInfo.singleStimCh;
    stimInfo=load('stimInfo.mat').stimInfo;
    extraInfo = load('extraInfo.mat').extraInfo;
    stimPulseArrs=load('stimPulseArrs.mat').arrs;
    
    fs = 30000;
    A = readmda('firings.curated.mda'); % read firings data
    primaryCh_vec = A(1,:); 
    timeStamp_samples = A(2,:); timeStamp_seconds = timeStamp_samples./fs;
    unit_vec = A(3,:);
    clear A; 
    
    if strcmp(extraInfo.stimType,'single')
        stimChannelInfo = channelInfo.singleStimCh;
        singleRowIdx= channelInfo.singleStimCh.singleRowIdx;
        single_stim_times = sorted_TS_current_singleCh(:,1);
    elseif strcmp(extraInfo.stimType,'paired')
        stimChannelInfo = channelInfo.pairedStimCh;
        pairedRowIdx= channelInfo.pairedStimCh.pairedRowIdx;
        paired_stim_times = sorted_TS_current_pairedCh(:,1);
    else
        stimChannelInfo = channelInfo.singleStimCh;
        singleRowIdx= channelInfo.singleStimCh.singleRowIdx;
        single_stim_times = sorted_TS_current_singleCh(:,1);
        paired_stim_times = sorted_TS_current_pairedCh(:,1);
    end
    
    clusterInfo = zeros(1,3);
    units = unique(unit_vec);
    nUNITS = numel(units);
    STIM_CHANNEL_NUM = numel(singleRowIdx);

    for i = 1:nUNITS
        unitIDX = units(i);
        cluster_info(i,1) = unitIDX; % not important 
        primaryCh = primaryCh_vec(find((unit_vec==unitIDX),1));
        cluster_info(i,2) = primaryCh; % map this to original channel name 
        cluster_info(i,3) = numel(find(cluster_info(:,2)==primaryCh));
    end

    filtered_data = readmda('filt.mda');

    allRowIdx=channelInfo.allCh.allRowIdx;
    all128OrderIdx=channelInfo.allCh.all128OrderIdx;
    allChLocations=channelInfo.allCh.allChLocations;
    allChNames=channelInfo.allCh.allChNames;
 
    h = figure(1); clf;
    h.Position = [61 194 1278 1084];
    
%     stimulating_ch = {};
    totalCh = numel(channelInfo.allCh.allRowIdx);
    stim_graph = zeros(totalCh,totalCh);

    for i = 3:nUNITS
        
        unit = units(i);
        primaryCh = cluster_info(i,2);
        IDforCh = cluster_info(i,3); % there can be more than 1 cluster/neuron for each primary channel 
        
        name = allChNames{primaryCh}; ch_idx = all128OrderIdx{primaryCh}; x_coord = allChLocations{primaryCh,1}; y_coord = allChLocations{primaryCh,2};

        figure(1); clf;
        t = sgtitle(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh)),'Interpreter','none');
        set(t, 'horizontalAlignment', 'right')

        % selectBoolean selects all indices in which the unit is equal to the unit of interest
        selectBoolean = (unit_vec == unit);
        primaryChannel = primaryCh_vec(find((selectBoolean==1),1));
        neuronTimestamps = timeStamp_seconds(selectBoolean)';
        neuronTimestampsSamples = timeStamp_samples (selectBoolean)';
        
        NUM_REPLICANTS = 10;
        currents = [0.5 1 2 5];
        nCU = numel(currents);
        % Keeps firings if timestamps fall in user-defined time boundaries
        [FIRING_DATA,post_stim_avg_window_1_array,post_stim_avg_window_2_array] = ...
                getEvents2(neuronTimestamps, stimPulseArrs.singleStimPulseArr, POSTSTIM_BLANK_TIME, POSTSTIM_AVG_WINDOW, ...
                PRE_STIM_WINDOW,POST_STIM_WINDOW,STIM_CHANNEL_NUM ,nCU,NUM_REPLICANTS,currents,all128OrderIdx);
            
%         [FIRING_DATA ,post_stim_avg_window_1_array,post_stim_avg_window_2_array] = ...
%             getEvents2(neuronTimestamps, stim_pulse_arr_paired, POSTSTIM_BLANK_TIME, POSTSTIM_AVG_WINDOW, ...
%             PRE_STIM_WINDOW,POST_STIM_WINDOW,nStimPairs,nCU,NUM_REPLICANTS,currents,ch_pair_idx);    
%     
        % plotting the raster plot on left-hand side
        % each row in the rasterplot is the spikes in a stimulation window
        stimChannelInfo = channelInfo.singleStimCh;
        generateRasterplot2(FIRING_DATA, stimChannelInfo, PRE_STIM_WINDOW,POST_STIM_WINDOW,STIM_CHANNEL_NUM,nCU,NUM_REPLICANTS,0);
        
        numPoints = (PRE_STIM_WINDOW + POST_STIM_WINDOW)/BIN_RESOLUTION ;
        time_hist = linspace(-PRE_STIM_WINDOW,POST_STIM_WINDOW,numPoints+1);
        subplot_tight(1,4,3)
        colors = get(gca,'colororder');
        
        for fig = 1:STIM_CHANNEL_NUM
        
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

                stimulated_ch = primaryCh; 
                stim_ch = singleRowIdx(fig);
                       
                if preStimMean*post_temp(1) > 2*preStimStd || preStimMean*post_temp(2) > 2*preStimStd
                    stim_graph(stim_ch,stimulated_ch) = stim_graph(stim_ch,stimulated_ch)+curr;
                else
                    stim_graph(stim_ch,stimulated_ch) = 0;
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

        %% plot the average waveform of "specific" bin in a small inset at the top right
        % the histogram outputs the histogram bin index FOR ALL UNITS WITHIN
        % THE LIMITS OF THE STIMULATION. Therefore, we need to pare down the
        % detected unit timestamps to just those within 5 msec of the
        % stimulation times.
    
        generateAvgPostStimWaveform(post_stim_avg_window_1_array,post_stim_avg_window_2_array,...
                             filtered_data,primaryCh,tBefore,tAfter)
        cd(image_save_path);
        export_fig(char(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh), '.png')),'-native')
        cd('..');

    end
    cd('..');
    
    G = digraph(stim_graph,allChNames);
    plot(G)
    title('ICMS-15 single stim digraph')
    
    
    
end
