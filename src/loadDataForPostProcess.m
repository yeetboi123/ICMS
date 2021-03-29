function [channelInfo,stimInfo,stimPulseArrs,primaryCh_vec,timeStamp_samples,unit_vec]     loadDataForPostProcess(processing_path)

    addpath(processing_path)
    
    channelInfo=load('channelInfo.mat');
    stimInfo=load('stimInfo.mat');
    stimPulseArrs=load('stimPulseArrs.mat');
    
    fs = 30000;
    A = readmda('firings.curated.mda'); % read firings data
    primaryCh_vec = A(1,:); 
    timeStamp_samples = A(2,:); timeStamp_seconds = timeStamp_samples./fs;
    unit_vec = A(3,:);
    clear A; 

    clusterInfo = zeros(1,3);
    units = unique(unit_vec);
    nUNITS = numel(units);
%     STIM_CHANNEL_NUM = numel(unique(sorted_TS_current_singleCh(:,3)));

    for i = 1:nUNITS
        unitIDX = units(i);
        cluster_info(i,1) = unitIDX; % not important 
        primaryCh = primaryCh_vec(find((unit_vec==unitIDX),1));
        cluster_info(i,2) = primaryCh; % map this to original channel name 
        cluster_info(i,3) = numel(find(cluster_info(:,2)==primaryCh_vec));
    end

    filtered_data = readmda('filt.mda');

%     single_stim_times = sorted_TS_current_singleCh(:,1);
%     paired_stim_times = sorted_TS_current_pairedCh(:,1);


end