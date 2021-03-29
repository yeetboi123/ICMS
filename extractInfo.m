function [channelInfo,varargout] = extractInfo(enabled_ch_names, new_order, stim_data, t, processing_path)
%extractInfo uses the outputs from loadData to produce intermediate files 

% extractInfo uses loaded data to create intermediate files that will be
% used for preprocessing the raw data 

% Write and save spatial map file for Mountainsort
coords = readmatrix('channel_maps/128channel4shank_location_data.csv');
coords_of_interest = coords(new_order,:);
writematrix(coords_of_interest, fullfile(processing_path ,'map.csv'));

% Extract stimulation pulse information 
[stim_info, stimType, stimmersPerPulse, changes, STIM_IDX_BEFORE_CHANGE] = extractStimInfo(stim_data);

% Write and save channel_info structure 
allRowIdx = (1:numel(new_order))';
channelInfo.allCh.allRowIdx = allRowIdx;
channelInfo.allCh.all128OrderIdx = enabled_ch_names;
channelInfo.allCh.allChNames = enabled_ch_names;
channelInfo.allCh.allChLocations = num2cell(coords_of_interest);

if strcmp(stimType,'single')
    singleRowIdx = unique(cell2mat(stim_info(1,:)));
    channelInfo.singleStimCh.singleRowIdx = singleRowIdx;
    channelInfo.singleStimCh.single128OrderIdx = new_order(singleRowIdx);
    channelInfo.singleStimCh.singleChNames = enabled_ch_names(singleRowIdx);
    channelInfo.singleStimCh.singleChLocations = num2cell(coords_of_interest(singleRowIdx,:));
%     channelInfo.singleStimCh = {singleRowIdx' single128OrderIdx' singleChNames singleChLocations'};
elseif strcmp(stimType,'paired')
    pairedRowIdx = unique(cell2mat(stim_info(1,:))','rows'); 
    channelInfo.pairedStimCh.pairedRowIdx = pairedRowIdx;
    channelInfo.pairedStimCh.paired128OrderIdx = new_order(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChNames = enabled_ch_names(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChLocations = makePairedStimCellArray(pairedRowIdx,coords_of_interest);
%     channelInfo.pairedStimCh = {pairedRowIdx' paired128OrderIdx' pairedChNames' pairedChLocations'};
else
    singleRowIdx = unique(cell2mat(stim_info(1,1:STIM_IDX_BEFORE_CHANGE)))';
    channelInfo.singleStimCh.singleRowIdx = singleRowIdx;  
    channelInfo.singleStimCh.single128OrderIdx = new_order(singleRowIdx)';
    channelInfo.singleStimCh.singleChNames = enabled_ch_names(singleRowIdx);
    channelInfo.singleStimCh.singleChLocations = num2cell(coords_of_interest(singleRowIdx,:))';
%     channelInfo.singleStimCh = {singleRowIdx' single128OrderIdx' singleChNames singleChLocations'};
    
    pairedRowIdx = unique(cell2mat(stim_info(1,STIM_IDX_BEFORE_CHANGE+1:end))','rows');
    channelInfo.pairedStimCh.pairedRowIdx = pairedRowIdx;
    channelInfo.pairedStimCh.paired128OrderIdx = new_order(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChNames = enabled_ch_names(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChLocations = makePairedStimCellArray(pairedRowIdx,coords_of_interest);
%     channelInfo.pairedStimCh = {pairedRowIdx paired128OrderIdx pairedChNames pairedChLocations};
end

if nargout > 1
    varargout = {stimType, stimmersPerPulse, changes, STIM_IDX_BEFORE_CHANGE};
end

end

function c = makePairedStimCellArray(pairedRowIdx,coords_of_interest)
    trials = size(pairedRowIdx,1);
    a = mat2cell(coords_of_interest(pairedRowIdx(:,1),:),ones(trials,1),2);
    b = mat2cell(coords_of_interest(pairedRowIdx(:,2),:),ones(trials,1),2);
    c = [a b];
end