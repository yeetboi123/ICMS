NUM_CORES = 12;
L = length(stim_times);
INTERVALS_PER_CORE = floor((L-1)/NUM_CORES);
new_data = zeros(size(raw,1),size(raw,2));
dsc_times = find(t == 0); dsc_times = dsc_times(2:end);
for CH_IND = 1:size(raw,1)
    tic
    CH_IND
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
            stim_idx{c} = t1:numel(stim_times);
        else
            data_arr{c} = sample(stim_times(t1):stim_times(t2)-1);
            stim_idx{c} = t1:t2;
        end    
    end


    yeet = [];
    parfor seg = 1:NUM_CORES
        ys = [];
        segment = data_arr{seg};
        sIDX = stim_idx{seg};
        norm_stim_times = stim_times(sIDX) - stim_times(sIDX(1)) + 1; 
        if seg == NUM_CORES
            sIDX(end) = [];
        end
        norm_dsc_times = dsc_times(sIDX) - stim_times(sIDX(1)) + 1;
        
        if seg == NUM_CORES
            IPC = numel(sIDX);
        else
            IPC = INTERVALS_PER_CORE;
        end

        for i = 1:IPC
            chunk1 = segment(norm_stim_times(i):norm_dsc_times(i)-1);
            chunk2 = segment(norm_dsc_times(i):norm_stim_times(i+1)-1);
            y1 = nt_detrend(chunk1',5,[],[],[],[],60);
            y2 = nt_detrend(chunk2',5,[],[],[],[],60);
            ys = [ys y1' y2'];
        end
        yeet = [yeet ys];
    end
    new_data(CH_IND,:) = [sample(1:stim_times(1)-1) yeet sample(stim_times(end):end)];
    toc
end

for i = 1:2
    subplot(2,1,i)
    plot(raw(i,1:end)')
    hold on;
    plot(new_data(i,1:end)')
    gridxy(stim_times)
    gridxy(dsc_times)
end


plot(raw(i,1:2e6)')
hold on;
plot(new_data(i,1:2e6)')
gridxy(stim_times)
gridxy(stim_times + 30)
% gridxy(dsc_times)

newer_data = new_data;
for ch = 1:size(filtered_data,1)
    ch_data = new_data(ch,:);
    for i = 1:numel(stim_times)
       
       newer_data(ch,stim_times(i):stim_times(i)+10) = zeros(1,11);
    end
end