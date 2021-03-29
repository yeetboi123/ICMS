function [raw, enabled_ch_names, new_order, stim_data, t, channelInfo, stimInfo, extraInfo]=preprocessData()
    addpath util\
    [raw, enabled_ch_names, new_order, stim_data, t]=loadData();


%     [enabled_ch_names new_order]=loadData;
    processing_path = uigetdir('Select directory to save intermediary files for preprocessing');
    
    [channelInfo,stimInfo,extraInfo] = extractInfo(enabled_ch_names, new_order, stim_data, t, processing_path);
%     extraInfo.stimType = {stimType, stimmersPerPulse, changes, STIM_IDX_BEFORE_CHANGE};
end