function [time_arr,current_arr] = getCathTimePaired(stim_data)

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

current_arr = {};
time_arr = {};
for i = 1:numel(cathodes)
    cath = cathodes(i);
    stim_times = find(stim_data(cath,:) < 0);
    phaseBegin = [1 (find(diff(stim_times) > 1) + 1)];
    phaseBegin = stim_times(phaseBegin);
    time_arr{i} = phaseBegin;
    current_arr{i} = -stim_data(cath,phaseBegin);
end

end