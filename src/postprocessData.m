function []=postprocessData(rerun,debug,filtered_data)
    %postprocessData(1,0,filtered_data)
    addpath(genpath('../src')) 
    addpath(genpath('../util'))
    
    processing_path = uigetdir('Select directory to read and write files for postprocessing'); 
    cd(processing_path)
    image_save_path = fullfile(processing_path,'output_images');
    if ~exist(image_save_path, 'dir')
        mkdir(image_save_path)
        addpath(image_save_path)
    end
    
    BIN_RESOLUTION = 0.01; % time resolution of histogram bin in seconds
    PRESTIM_BLANK_TIME= 0.000 ; % zero out output t seconds before stim
    POSTSTIM_BLANK_TIME = 0.000; % zero out output t seconds after stim
    POSTSTIM_AVG_WINDOW = 0.005; % window of post-stim average waveform 
    PRE_STIM_WINDOW = 0.4; % window before and after stim for plotting rasterplot 
    POST_STIM_WINDOW = 0.4; % window before and after stim for plotting rasterplot 
    currents = [0.5 1 2 5]; % IMPORTANT: stimulation currents used!
    
    if ~rerun
        filtered_data = readmda('filt.mda'); % read filtered data
        moveToBaseWorkspace(filtered_data);
    end
    
    A = readmda('firings.curated.mda'); % read curated firings data from Mountainsort
    % load intermediate files
    channelInfo=load('channelInfo.mat').channelInfo;
    stimChannelInfo = channelInfo.singleStimCh;
    stimInfo=load('stimInfo.mat').stimInfo;
    extraInfo = load('extraInfo.mat').extraInfo;
    stimPulseArrs=load('stimPulseArrs.mat');
    coordinates = channelInfo.allCh.allChLocations;
    
    if strcmp(extraInfo.stimType,'single')
        stimChannelInfo = channelInfo.singleStimCh;
        singleRowIdx= channelInfo.singleStimCh.singleRowIdx;
        stim_pulse_arr = stimPulseArrs.singleStimPulseArr;
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
    
    % info of all channels recorded 
    allRowIdx=channelInfo.allCh.allRowIdx;
    all128OrderIdx=channelInfo.allCh.all128OrderIdx;
    allChLocations=channelInfo.allCh.allChLocations;
    allChNames=channelInfo.allCh.allChNames;
    
    % assign rows from firings to variables
    fs = 30000; 
    primaryCh_vec = A(1,:); % primary channel of unit (channel with highest amplitude of unit), 128-order indexing
    timeStamp_samples = A(2,:); timeStamp_seconds = timeStamp_samples./fs;
    unit_vec = A(3,:);
    clear A; 

    units = unique(unit_vec); % indices of units with respect to Mountainsort cluster indices
    nUNITS = numel(units); % number of units
    nSCH = numel(singleRowIdx); % number of stim channels
    
    cluster_info = zeros(nUNITS,3);
    for i = 1:nUNITS
        unitIDX = units(i);
        cluster_info(i,1) = unitIDX; % not important 
        primaryCh = primaryCh_vec(find((unit_vec==unitIDX),1));
        cluster_info(i,2) = primaryCh; % map this to original channel name 
        cluster_info(i,3) = numel(find(cluster_info(:,2)==primaryCh));
    end

    
    for i = 1:numel(units)

        totalCh = numel(allRowIdx);
        stim_graph = zeros(totalCh,totalCh); % create empty stim graph 
        
        unit = units(i);
        primaryCh = cluster_info(i,2);
        IDforCh = cluster_info(i,3); % there can be more than 1 cluster/neuron for each primary channel 
        
        name = allChNames{primaryCh}; x_coord = allChLocations{primaryCh,1}; y_coord = allChLocations{primaryCh,2};

        figure('units','normalized','outerposition',[0 0 1 1]); clf;
        t = sgtitle(strcat('Unit',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh)),'Interpreter','none');
        set(t, 'horizontalAlignment', 'right')
        
        % selectBoolean selects all indices in which the unit is equal to the unit of interest
        selectBoolean = (unit_vec == unit);
        neuronTimestamps = timeStamp_seconds(selectBoolean)';
        neuronTimestampsSamples = timeStamp_samples (selectBoolean)';
        
        NUM_REPLICANTS = size(stim_pulse_arr{1,1},2);  
        nCU = numel(currents);
        % Keeps firings if timestamps fall in user-defined time boundaries
        [FIRING_DATA,post_stim_avg_window_1_array,post_stim_avg_window_2_array] = ...
                getEvents2(neuronTimestamps, stimPulseArrs.singleStimPulseArr, POSTSTIM_BLANK_TIME, POSTSTIM_AVG_WINDOW, ...
                PRE_STIM_WINDOW,POST_STIM_WINDOW,nSCH,nCU,NUM_REPLICANTS,currents,all128OrderIdx);
            
%         [FIRING_DATA ,post_stim_avg_window_1_array,post_stim_avg_window_2_array] = ...
%             getEvents2(neuronTimestamps, stim_pulse_arr_paired, POSTSTIM_BLANK_TIME, POSTSTIM_AVG_WINDOW, ...
%             PRE_STIM_WINDOW,POST_STIM_WINDOW,nStimPairs,nCU,NUM_REPLICANTS,currents,ch_pair_idx);    
%     
        % plotting the raster plot on left-hand side
        % each row in the rasterplot is the spikes in a stimulation window
        stimChannelInfo = channelInfo.singleStimCh;
        generateRasterplot2(FIRING_DATA, stimChannelInfo,PRE_STIM_WINDOW,POST_STIM_WINDOW,nSCH,nCU,NUM_REPLICANTS,0);
        
        numPoints = (PRE_STIM_WINDOW + POST_STIM_WINDOW)/BIN_RESOLUTION ;
        time_hist = linspace(-PRE_STIM_WINDOW,POST_STIM_WINDOW,numPoints+1);
        subplot_tight(1,4,3,[0.05,0.03]) % make sure margins same as those of rasterplot subplot for alignment
        colors = get(gca,'colororder');
        
        vert_space = nCU * NUM_REPLICANTS;
        
        
        for j = 1:nSCH % iterate through stimulating channels
            
            height_offset = vert_space * (j-1) + 20; % shift PSTH every stimulating channel 
            combinedPSTH = [];
            
            curr_firings = [];
            activating_currents = [];
            for curr = 1:nCU
                % take average of replicants
                for rep = 1:NUM_REPLICANTS
                    temp = [FIRING_DATA{j,curr}{rep,1}{1,2}]; % check this
                    curr_firings = [curr_firings; temp];
                end    
               [allCounts,~,~] = histcounts(curr_firings,time_hist);
               [preStim_count, postStim_count, counts,preStimMean,preStimStd] = getHistLimits(numPoints, allCounts,time_hist);

                postStim_limit = floor(0.02/BIN_RESOLUTION); % 20 ms window after stim considered post stim period 
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
                stim_ch = singleRowIdx(j);

                activations = find(postStim > (preStimMean+3*preStimStd)/preStimMean); % bins   
                
                if ~isempty(activations)
                    activating_currents = [activating_currents curr];
                    actIDXs = activations(1); % consecutive activation bin array 
                    idx = 2;
                    while 1
                        if (idx > length(activations))
                            break
                        % store if activations next to each other
                        elseif (activations(idx) - activations(idx-1) == 1)
                            actIDXs = [actIDXs activations(idx)]; % append
                            idx = idx + 1;
                        else
                            break
                        end
                    end
                    actIDXs = actIDXs + postStim_count - 1;
                end
                
                
                % overlay a green bar over the activations
                
%                 x_offset = preStim_count + ((length(time_hist)-1) * curr);
%                 x1 = x_offset + time_hist(actInds(1)); y1 = 
                
%                 v_j = [time_hist(actInds(1)) 0; time_hist(actInds(end)+1) 0; ...
%                     time_hist(actInds(end)+1) 1.25*max(counts); time_hist(actInds(1)) 1.25*max(counts)];
                
%                 if preStimMean*post_temp(1) > 2*preStimStd || preStimMean*post_temp(2) > 2*preStimStd
%                     stim_graph(stim_ch,stimulated_ch) = stim_graph(stim_ch,stimulated_ch)+curr;
%                     activating_currents = [activating_currents curr];
%                 else
%                     stim_graph(stim_ch,stimulated_ch) = 0;
%                 end

            end
%             height = vert_space * (j-1) + 20;
            plot((30*combinedPSTH)+120 * (j-1),'LineWidth',1,'color',colors(7,:))
%             patch('Faces',f,'Vertices',v,'FaceColor','green','FaceAlpha',.3,'EdgeColor','none'); % f will be repeating array of 1234
            if ~isempty(activating_currents)
                min_act_curr = min(activating_currents);
                stim_graph(stim_ch,stimulated_ch) = stim_graph(stim_ch,stimulated_ch) + ...
                                                    min_act_curr; % store minimal activating current in graph 
            end
            
            hold on;
            gridxy(preStim_count + ((length(time_hist)-1) * (0:3)),'LineStyle',':','linewidth',0.25);
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
        ax1=generateAverageWaveform(neuronTimestampsSamples,...
            filtered_data,primaryCh,tBefore,tAfter);

        %% plot the average waveform of "specific" bin in a small inset at the top right
        % the histogram outputs the histogram bin index FOR ALL UNITS WITHIN
        % THE LIMITS OF THE STIMULATION. Therefore, we need to pare down the
        % detected unit timestamps to just those within 5 msec of the
        % stimulation times.
    
        [ax2,ax3]=generateAvgPostStimWaveform(post_stim_avg_window_1_array,post_stim_avg_window_2_array,...
                             filtered_data,primaryCh,tBefore,tAfter);
        
        generateActivationDigraph(stim_graph,stimulated_ch,nCU,coordinates)
        
        linkaxes([ax1 ax2 ax3],'xy')  
        shift = 0.008;
        pos1 = get(ax1, 'Position'); pos2 = get(ax2, 'Position'); pos3 = get(ax3, 'Position'); 
        posnew1 = pos1; posnew1(2) = posnew1(2) + shift;
        posnew2 = pos2; posnew2(2) = posnew2(2) + shift; 
        posnew3 = pos3; posnew3(2) = posnew3(2) + shift; 
        set(ax1, 'Position', posnew1); set(ax2, 'Position', posnew2); set(ax3, 'Position', posnew3)
        
        cd(image_save_path);
        if ~debug
            export_fig(char(strcat('Cluster',{'  '},name,{'  '},num2str(x_coord) ,',', num2str(y_coord), '  #', num2str(IDforCh), '.png'))...
            ,'-nocrop','-p0.03','-native')
        end
        cd('..');

    end
    
    generateStimulationDigraph(stim_graph,nSCH,nCU,singleRowIdx,coordinates,stimChannelInfo,processing_path)
    
    
    cd('..');
    close all;
       
    
end

% function moveToBaseWorkspace(variable)
% % move_to_base_workspace(variable)
% %
% % Move variable from function workspace to base MATLAB workspace so
% % user will have access to it after the program ends.
% variable_name = inputname(1);
% assignin('base', variable_name, variable);
% 
% end