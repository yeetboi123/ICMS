function stimPulseArrs = makeStimPulseArrays(channelInfo,stimInfo,extraInfo,processing_path)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
stimIdxBeforeChange = extraInfo.stimIdxBeforeChange;
stimType = extraInfo.stimType;

if strcmp(stimType,'single')
    single_stim_info = stimInfo;
elseif strcmp(stimType,'paired')
    paired_stim_info = stimInfo;
else
    % split the stimInfo array into single stimulation and paired
    % stimulation arrays 
    single_stim_info = cell2mat(stimInfo(:,1:stimIdxBeforeChange));   
    paired_stim_info = stimInfo(:,stimIdxBeforeChange+1:end);
end

if strcmp(stimType,'single') || strcmp(stimType,'combined')
    
    single_stim_row_idx = channelInfo.singleStimCh.singleRowIdx; 
    single_stim_chs = single_stim_info(1,:);
    single_stim_times = single_stim_info(2,:);
    single_stim_currents = single_stim_info(3,:);
    currents = flip(unique(single_stim_currents));
    
    getCurrentInd = @(x) find(abs(single_stim_currents - x) < 0.001);
    getChInd = @(x) find(single_stim_chs == x);
    getStimPulse = @(ch,i) single_stim_times(intersect(getCurrentInd(i),getChInd(ch)));
%     NUM_REPLICANTS = numel(getStimPulse(single_stim_idx(1),currents(1)));
    
    nCH = numel(single_stim_row_idx); nCU = numel(currents); 
    singleStimPulseArr = cell(nCH,nCU);
    for i = 1:nCH
        for j = 1:nCU
                stim_ch = single_stim_row_idx(i); stim_current = currents(j);
                singleStimPulseArr{i,j} = getStimPulse(stim_ch,stim_current);
                
% Unit test
%                 for k = 1:length(singleStimPulseArr{i,j})
%                     ts = singleStimPulseArr{i,j}(k);
%                     a = find(stim_data(:,ts) < 0);
%                     assert(a == single_stim_row_idx(i))
%                     assert(stim_data(a,ts) == currents(j))
%                 end                
        end
    end
    stimPulseArrs.singleStimPulseArr = singleStimPulseArr;
end

% ts = singleStimPulseArr{2,2}
% plot(stim_data(single_stim_row_idx(2),ts(1)-100:ts(1)+100))

if strcmp(stimType,'paired') || strcmp(stimType,'combined')
    
    paired_stim_row_idx = channelInfo.pairedStimCh.pairedRowIdx; 
    paired_stim_chs = cell2mat(paired_stim_info(1,:));
    paired_stim_times = cell2mat(paired_stim_info(2,:));
    paired_stim_currents = cell2mat(paired_stim_info(3,:));
    currents = flip(unique(paired_stim_currents));
    
    getCurrentInd = @(x) find(abs(paired_stim_currents - x) < 0.001);
    getChInd = @(x) find(ismember(paired_stim_chs',x,'rows'));
    getStimPulse = @(ch,i) paired_stim_times(intersect(getCurrentInd(i),getChInd(ch)));

    nCH = size(paired_stim_row_idx,1); nCU = numel(currents); 
    
    pairedStimPulseArr = cell(nCH,nCU);
    for i = 1:nCH
        for j = 1:nCU
                stim_ch = paired_stim_row_idx(i,:); stim_current = currents(j);
                pairedStimPulseArr{i,j} = getStimPulse(stim_ch,stim_current);
        end
    end
    stimPulseArrs.pairedStimPulseArr = pairedStimPulseArr;
end

% ts = pairedStimPulseArr{5,4}
% plot(stim_data(paired_stim_row_idx(5,:),ts(1)-100:ts(1)+100)')
save(fullfile(processing_path,'stimPulseArrs.mat'),'-struct','stimPulseArrs');

end
