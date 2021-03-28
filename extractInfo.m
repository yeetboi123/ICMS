function extractInfo(raw, enabled_ch_names, new_order, stim_data, t, processing_path, SINGLE_AND_OR_PAIRED)
%extractInfo uses the outputs from loadData to produce intermediary files 

% extractInfo uses loaded data to create intermediary files that will be
% used for preprocessing the raw data 

% Write and save spatial map file for Mountainsort
coords = readmatrix('channel_maps/128channel4shank_location_data.csv');
coords_of_interest = coords(new_order,:);
writematrix(coords_of_interest, fullfile(processing_path ,'map.csv'));

% Find stimulation channels for single stimulation and/or paired stimulation 
switch SINGLE_AND_OR_PAIRED
    case 1
        disp('Generating intermediary files for single channel stimulation data...')
        single_stim_ch = find(any(stim_data,2) == 1);
    case 2 
        disp('Generating intermediary files for paired channel stimulation data...')
        paired_stim_ch = find(any(stim_data,2) == 1);
    case 3
        disp('Generating intermediary files for single and paired channel stimulation data...')
        stim_ch = find(any(stim_data,2) == 1);
end

% Write and save channel_info structure 
channelInfo.allCh = {};
channelInfo.recordingCh = {};
channelInfo.singleStimCh = {};
channelInfo.pairedStimCh = {};

allRowIdx = 1:numel(new_order)';
all128OrderIdx = enabled_ch_names';
allChNames = enabled_ch_names';
allChLocations = num2cell(coords_of_interest)';

channelInfo.allCh = {allRowIdx all128OrderIdx allChNames allChLocations};
channelInfo.singleStimCh = find(

end

function [time_arr,current_arr,stim_order] = getCathTimeNew(stim_data)

% the stimulation channels are those with any value that is not zero
stim_chs=find(any(stim_data,2) == 1);

% here we separate the cathodes (negative phase first) and anodes (positive
% phase first). Timing will be set to the cathodal pulse
cathodes = []; anodes = [];
for i = 1:numel(stim_chs)
    firstStim = find((stim_data(stim_chs(i),:))~=0,1);
    firstStim = stim_data(stim_chs(i),firstStim);
    if (firstStim > 0)
        anodes = [anodes, stim_chs(i)];
    elseif (firstStim < 0)
        cathodes = [cathodes, stim_chs(i)];
    end
end

current_arr = [];
time_arr = [];
for i = cathodes
    stim_times = find(stim_data(i,:) < 0);
    phaseBegin = [1 (find(diff(stim_times) > 1) + 1)];
    phaseBegin = stim_times(phaseBegin);
    time_arr = [time_arr; phaseBegin];
    current_arr = [current_arr; -stim_data(i,phaseBegin)];
end

stim_order = [];
for j = 1:numel(stim_times)
    stim_time = stim_times(j);
    stim_ch = find(stim_data(:,stim_time) < 0);
    stim_order = [stim_order stim_ch];
end