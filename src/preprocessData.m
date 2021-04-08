function [raw, enabled_ch_names, new_order, stim_data, t, channelInfo, stimInfo, extraInfo]=preprocessData()
    addpath(genpath('src')) 
    addpath(genpath('util'))
    processing_path = uigetdir('Select directory to save intermediary files for preprocessing');
    addpath(processing_path); cd(processing_path)
%     [raw, enabled_ch_names, new_order, stim_data, t] = loadData(); %#ok<ASGLU> % output all for new data 
    [enabled_ch_names, new_order, stim_data, t] = loadData(); % load without raw data for debugging purposes
    cd ..
    [channelInfo,stimInfo,extraInfo] = extractInfo(enabled_ch_names, new_order, stim_data, t, processing_path);
    stimPulseArrs = makeStimPulseArrays(channelInfo,stimInfo,extraInfo,processing_path);
    
    % preprocess raw data
    replicant_array = createReplicantArray(raw,stimInfo,channelInfo,stimPulseArrs,300,0);
%     plot(output_data(1,1:1e6)); hold on; plot(raw(1,1:1e6))
    [~,output_data] = artifactRemoveData(replicant_array,raw,stimPulseArrs);
    
    
    detrended = detrendData2(output_data,stimInfo,t,nSamples,processing_path);
    plot(detrended(4,1:1e6)); hold on; plot(raw(4,1:1e6))
    writemda32(detrended,'preprocessed_data.mda');
end

