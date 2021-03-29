function [raw, enabled_ch_names, new_order, stim_data, t, channelInfo, stimInfo, extraInfo]=preprocessData()
    addpath util\
    [raw, enabled_ch_names, new_order, stim_data, t]=loadData();
    processing_path = uigetdir('Select directory to save intermediary files for preprocessing');
    [channelInfo,stimInfo,extraInfo] = extractInfo(raw, enabled_ch_names, new_order, stim_data, t, processing_path);
    arrs = makeStimPulseArrays(channelInfo,stimInfo,extraInfo,processing_path);
    
    % preprocess raw data
    detrended = detrendData(raw-mean(raw),stimInfo,t,processing_path);
end

