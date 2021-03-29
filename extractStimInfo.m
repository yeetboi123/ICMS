function [stim_info, stimType, varargout] = extractStimInfo(stim_data)

    [~,t]=find(stim_data ~= 0);
    tt=unique(t); cth_phz = tt(1:6:end);

    stim_info = cell(3,numel(cth_phz));
    for k = 1:numel(cth_phz)
        stim_time = cth_phz(k);
        stim_info{1,k} = find(stim_data(:,stim_time) < 0); % stim channel(s) 
        stim_info{2,k} = cth_phz(k); % beginning of cathodic phase
        stimmers = stim_info{1,k}; 
        stim_info{3,k} = stim_data(stimmers(1),cth_phz(k)); % stim current
    end
    
    stimmersPerPulse = cellfun(@(x) numel(x),stim_info(1,:));
    changes = nnz(diff(stimmersPerPulse));
    STIM_IDX_BEFORE_CHANGE = find(diff(stimmersPerPulse)==1);

    switch changes
        case 0
            if stimmersPerPulse(1) == 1
                stimType = 'single';
                disp('Single stimulation only file detected!');
            elseif stimmersPerPulse(1) == 2
                stimType = 'paired';
                disp('Paired stimulation only file detected!');
            else
                stimType = 'error';
                disp('Error! More than 2 stimulation pulses in trials! Faulty data acquisition!')
                return
            end
        case 1
            if stimmersPerPulse(1) == 1
                stimType = 'combined';
                disp('Proper combined stimulation file detected!');
            elseif stimmersPerPulse(1) == 2
                stimType = 'error';
                disp('Error! Paired stim trials before single stim trials!');
                return
            end
        otherwise
            stimType = 'error';
            disp('Error! Mixed stim trial types detected! Faulty data acquisition!');
            return
    end
          
    if nargout > 2
        varargout = {stimmersPerPulse changes STIM_IDX_BEFORE_CHANGE};
    end
    
end
