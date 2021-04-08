function [stim_info, extraInfo] = extractStimInfo(stim_data,t)

    [~,col]=find(stim_data ~= 0);
    tt=unique(col); 
    if tt(1) == 1 % stim pulse at sample 1 can occur
        cth_phz = tt(2:6:end);
    else
        cth_phz = tt(1:6:end);
    end
    
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
    STIM_IDX_BEFORE_CHANGE = find(abs(diff(stimmersPerPulse))==1);
    
    % check to make sure number of trials == number of stim pulses
    try 
        TRIAL_NUM = numel(find(t == 0)); 
        STIM_NUM = numel(cth_phz);
        assert(TRIAL_NUM == STIM_NUM);
    catch 
        error('Different number of trials (%d) and stimulation pulses (%d)! Faulty data acquisition!',...
            TRIAL_NUM,STIM_NUM')
    end
    
    switch changes
        case 0
            if stimmersPerPulse(1) == 1
                stimType = 'single';
                disp('Single stimulation only file detected!');
                [counts,stimch] = groupcounts([stim_info{1,1:end}]');
                stim_info = cell2mat(stim_info);
                if all(counts == counts(1))
                    disp('Identical number of pulses per single stim channel!')
                    fprintf('%d pulses per stim channel!\n', counts(1))
                else
                    fprintf('%d pulses for stim channel %d!\n', [counts;stimch])
                    fprintf('\n')
                    error('Different number of pulses per single stim channel! Faulty data acquisition!')
                end
            elseif stimmersPerPulse(1) == 2
                stimType = 'paired';
                disp('Paired stimulation only file detected!');
            else
                error('Error! More than 2 stimulation pulses in trials! Faulty data acquisition!');
            end
        case 1
            if stimmersPerPulse(1) == 1
                stimType = 'combined';
                disp('Proper combined stimulation file detected!');
                [counts,stimch] = groupcounts([stim_info{1,1:STIM_IDX_BEFORE_CHANGE}]');
                if all(counts == counts(1))
                    disp('Identical number of pulses per single stim channel!');
                    fprintf('%d pulses per stim channel!\n', counts(1))
                else
                    fprintf('%d pulses for stim channel %d!\n', [counts;stimch])
                    fprintf('\n')
                    error('Different number of pulses per single stim channel! Faulty data acquisition!');   
                end
            elseif stimmersPerPulse(1) == 2
                stimType = 'error';
                warning('Warning! Paired stim trials before single stim trials!');
                [counts,stimch] = groupcounts([stim_info{1,STIM_IDX_BEFORE_CHANGE+1:end}]');
                if all(counts == counts(1))
                    disp('Identical number of pulses per single stim channel!')
                    fprintf('%d pulses per stim channel!\n', counts(1))
                else
                    warning('Different number of pulses per single stim channel! Faulty data acquisition!')                
                    fprintf('%d pulses for stim channel %d!\n', [counts;stimch])
                    fprintf('\n')
                end
                    
            end
        otherwise
            error('Error! Mixed stim trial types detected! Faulty data acquisition!');
    end
     
    extraInfo.stimType = stimType;
    extraInfo.stimChPerPulse = stimmersPerPulse;
    extraInfo.changesInStimChPerPulse = changes;
    extraInfo.stimIdxBeforeChange = STIM_IDX_BEFORE_CHANGE;
    
end
