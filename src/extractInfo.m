function [channelInfo,stimInfo,extraInfo] = extractInfo(raw, enabled_ch_names, new_order, stim_data, t, processing_path)
%extractInfo uses the outputs from loadData to produce intermediate files 

% extractInfo uses loaded data to create intermediate files that will be
% used for preprocessing the raw data 

% Write and save spatial map file for Mountainsort
coords = readmatrix('channel_maps/128channel4shank_location_data.csv');
coords_of_interest = coords(new_order,:);
writematrix(coords_of_interest, fullfile(processing_path ,'map.csv'));

% Extract stimulation pulse information 
[stimInfo, extraInfo] = extractStimInfo(stim_data,t);
stimType = extraInfo.stimType;
STIM_IDX_BEFORE_CHANGE = extraInfo.stimIdxBeforeChange;

if strcmp(stimType,'error') 
    channelInfo = [];
    return
end
    

% Write and save channel_info structure 
allRowIdx = (1:numel(new_order))';
channelInfo.allCh.allRowIdx = allRowIdx;
channelInfo.allCh.all128OrderIdx = new_order(allRowIdx);
channelInfo.allCh.allChNames = enabled_ch_names;
channelInfo.allCh.allChLocations = num2cell(coords_of_interest);

if strcmp(stimType,'single')
    singleRowIdx = unique(cell2mat(stimInfo(1,:)));
    channelInfo.singleStimCh.singleRowIdx = singleRowIdx;
    channelInfo.singleStimCh.single128OrderIdx = new_order(singleRowIdx);
    channelInfo.singleStimCh.singleChNames = enabled_ch_names(singleRowIdx);
    channelInfo.singleStimCh.singleChLocations = num2cell(coords_of_interest(singleRowIdx,:));
%     channelInfo.singleStimCh = {singleRowIdx' single128OrderIdx' singleChNames singleChLocations'};
elseif strcmp(stimType,'paired')
    pairedRowIdx = unique(cell2mat(stimInfo(1,:))','rows'); 
    channelInfo.pairedStimCh.pairedRowIdx = pairedRowIdx;
    channelInfo.pairedStimCh.paired128OrderIdx = new_order(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChNames = enabled_ch_names(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChLocations = makePairedStimCellArray(pairedRowIdx,coords_of_interest);
%     channelInfo.pairedStimCh = {pairedRowIdx' paired128OrderIdx' pairedChNames' pairedChLocations'};
else
    singleRowIdx = unique(cell2mat(stimInfo(1,1:STIM_IDX_BEFORE_CHANGE)))';
    channelInfo.singleStimCh.singleRowIdx = singleRowIdx;  
    channelInfo.singleStimCh.single128OrderIdx = new_order(singleRowIdx)';
    channelInfo.singleStimCh.singleChNames = enabled_ch_names(singleRowIdx);
    channelInfo.singleStimCh.singleChLocations = num2cell(coords_of_interest(singleRowIdx,:));
%     channelInfo.singleStimCh = {singleRowIdx' single128OrderIdx' singleChNames singleChLocations'};
    
    pairedRowIdx = unique(cell2mat(stimInfo(1,STIM_IDX_BEFORE_CHANGE+1:end))','rows');
    channelInfo.pairedStimCh.pairedRowIdx = pairedRowIdx;
    channelInfo.pairedStimCh.paired128OrderIdx = new_order(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChNames = enabled_ch_names(pairedRowIdx);
    channelInfo.pairedStimCh.pairedChLocations = makePairedStimCellArray(pairedRowIdx,coords_of_interest);

end

checkpoint_data.raw = raw;
checkpoint_data.enabled_ch_names = enabled_ch_names;
checkpoint_data.new_order = new_order;
checkpoint_data.stim_data = stim_data;
checkpoint_data.t = t;
checkpoint_data.processing_path = processing_path;
checkpoint_data.channelInfo = channelInfo;
checkpoint_data.stimInfo = stimInfo;
checkpoint_data.extraInfo = extraInfo;
move_to_base_workspace(checkpoint_data); 

save(fullfile(processing_path,'channelInfo.mat'),'channelInfo');
save(fullfile(processing_path,'stimInfo.mat'),'stimInfo');
save(fullfile(processing_path,'extraInfo.mat'),'extraInfo');

end

function c = makePairedStimCellArray(pairedRowIdx,coords_of_interest)
    trials = size(pairedRowIdx,1);
    a = mat2cell(coords_of_interest(pairedRowIdx(:,1),:),ones(trials,1),2);
    b = mat2cell(coords_of_interest(pairedRowIdx(:,2),:),ones(trials,1),2);
    c = [a b];
end

function move_to_base_workspace(variable)
% move_to_base_workspace(variable)
%
% Move variable from function workspace to base MATLAB workspace so
% user will have access to it after the program ends.
variable_name = inputname(1);
assignin('base', variable_name, variable);

end