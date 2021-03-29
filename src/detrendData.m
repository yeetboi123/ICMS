function detrended = detrendData(raw,stimInfo,t,processing_path)
    addpath util\NoiseTools 
    addpath util\mdaio
    
    c = parcluster('local');
    NUM_CORES = c.NumWorkers;
    stim_times = cell2mat(stimInfo(2,:));
    
    L = numel(stim_times);
    INTERVALS_PER_CORE = max([floor(L/NUM_CORES) 1]); % at least 1 interval/core
    REMAINING_INTERVALS = L - (NUM_CORES - 1) * INTERVALS_PER_CORE;
    interval_assignments = zeros(NUM_CORES,1);
    for c = 1:NUM_CORES
        if c ~= NUM_CORES
            interval_assignments(c) = INTERVALS_PER_CORE;
        else
            interval_assignments(c) = REMAINING_INTERVALS;
        end
    end
    assert(sum(interval_assignments) == L,...
        'Interval assignments (%d) not adding up to number of stim pulses! (%d)',...
        sum(interval_assignments), L)
    
    dsc_times = find(t == 0); dsc_times = dsc_times(2:end);    
    detrended = zeros(size(raw,1),size(raw,2));

    for CH_IND = 1:size(raw,1)
        tic
        fprintf('Detrending channel %d...',CH_IND)
        sample = raw(CH_IND,1:end);

        % break up data into segments for parallel processing
        data_arr = cell(NUM_CORES,1);
        stim_idx = cell(NUM_CORES,1);
        prev_t2 = 1;
        for c = 1:NUM_CORES
            t1 = prev_t2;
            t2 = c * INTERVALS_PER_CORE + 1;
            prev_t2 = t2;  
            if c == NUM_CORES 
                data_arr{c} = sample(stim_times(t1):stim_times(end)-1);
                data_arr{c} = sample(stim_times(t1):end); % last segment stretches to the end
                stim_idx{c} = t1:numel(stim_times);
            else
                data_arr{c} = sample(stim_times(t1):stim_times(t2)-1);
                stim_idx{c} = t1:t2;
            end    
        end

        new_sample = cell(1,NUM_CORES);
        parfor seg = 1:NUM_CORES
            ys = [];
            segment = data_arr{seg};
            sIDX = stim_idx{seg};
            norm_stim_times = stim_times(sIDX) - stim_times(sIDX(1)) + 1; % normalize

            if seg == NUM_CORES
                norm_dsc_times = dsc_times(sIDX(1:numel(sIDX)-1)) - stim_times(sIDX(1)) + 1;
                IPC = numel(sIDX);
            else
                norm_dsc_times = dsc_times(sIDX) - stim_times(sIDX(1)) + 1;
                IPC = INTERVALS_PER_CORE;
            end

            for i = 1:IPC
                if seg == NUM_CORES && i == IPC
                    chunk1 = segment(norm_stim_times(i)+30:end);
                    y1 = beginBlank(chunk1);
                    ys = [ys y1'];
                else
                    chunk1 = segment(norm_stim_times(i)+30:norm_dsc_times(i)-1); 
                    chunk2 = segment(norm_dsc_times(i):norm_stim_times(i+1)-1);
                    y1 = beginBlank(chunk1);
                    y2 = nt_detrend(chunk2',5,[],[],[],[],60);
                    y2(1:15) = ones(1,15);
                    ys = [ys y1' y2'];
                end
            end
            
            new_sample{seg} = ys;
        end
        combined = cell2mat(new_sample);
        
%         plot(detrended(48,1e6:3e6));
%         hold on;
%         plot(raw(48,1e6:3e6));

        detrended(CH_IND,:) = [zeros(1,stim_times(1)-1) combined];
        toc
    end
    
    writemda32(detrended,fullfile(processing_path,'preprocessed_data.mda'));
    
end



