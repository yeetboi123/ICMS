% Artifact mitigation script

%% 1. Import raw data and 
addpath mdaio

[raw,stim_data,new_order,save_path] = createMappingNEW;
clear amp_settle_data charge_recovery_data compliance_limit_data t 

nCH = size(raw,1);
nSamples = size(raw,2);
STIM_PER_BLOCK = 1;
NUM_REPLICANTS = 11;
nCU = 4;
currents = [0.5 1 2 5];

[time_arr,current_arr] = getCathTimeNew(stim_data); time_arr_s = time_arr/3e4;
time_vector_unsorted = reshape(time_arr', 1, []); current_vector = reshape(current_arr', 1, []);
temp = repmat(new_order',1,nCU*NUM_REPLICANTS)'; channel_vector = temp(:)';

[B,current] = sort(time_vector_unsorted','ascend');
sorted_TS_current = [B current_vector(current)' channel_vector(current)'];
save(fullfile(save_path,'sorted_TS_current.mat'),'sorted_TS_current')
 %%  2. Channel mean subtraction (for noisy data)
% mean_sg = mean(raw);
% med_sg = median(raw);
% raw = raw - mean_sg; clear mean_sg
% t = (12e4:24e4)/3e1;
% subplot 211
% plot(t,raw(8,12e4:24e4)); xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Raw data')
% subplot 212
% plot(t,raw(8,12e4:24e4)- mean_sg(12e4:24e4)); xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Mean-subtracted data')
% sgtitle('ICMS15: Raw data vs. mean-subtracted data')
% 
% subplot 211
% plot(raw(:,1:30e4)')
% subplot 212
% plot((raw(:,1:30e4) - med_sg(1:30e4))')

%% 2. Remove strange discontinuities seen in Intan recordings 
load('sorted_TS_current.mat');
stim_times = sort(sorted_TS_current(:,1));
dsc_times = []; % array containing discontinuity timestamps

data = raw(1,:);
for stim_idx = 1:numel(stim_times)-1
    start_idx = stim_times(stim_idx) + 100;   
    end_idx = stim_times(stim_idx + 1) - 100;
    y = abs(diff(data(start_idx:end_idx)));
    [~,ind] = max(y);
    corrected_ind = ind + start_idx;
    dsc_times = [dsc_times corrected_ind];
end

%% 3. Detrend (takes a while)
addpath NoiseTools\

NUM_CORES = 12;
L = length(stim_times);
INTERVALS_PER_CORE = floor((L-1)/NUM_CORES);
detrended = zeros(size(raw,1),size(raw,2));

% for CH_IND = 1:size(raw,1)
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
            chunk1 = segment(norm_stim_times(i):norm_dsc_times(i)-1); % ignore first millisecond after stim
            chunk2 = segment(norm_dsc_times(i):norm_stim_times(i+1)-1);
            y1 = nt_detrend(chunk1',5,[],[],[],[],60);
            y2 = nt_detrend(chunk2',5,[],[],[],[],60);
            ys = [ys y1' y2']; % replace first ms after stim with zeros
        end
        yeet = [yeet ys];
    end
    detrended(CH_IND,:) = [sample(1:stim_times(1)-1) yeet sample(stim_times(end):end)];
    toc
end

load('sorted_TS_current.mat');
stim_times = sort(sorted_TS_current(:,1));

blanked = detrended;
for ch = 1:size(detrended,1)
    ch_data = detrended(ch,:);
    for i = 1:numel(stim_times)
       blanked(ch,stim_times(i):stim_times(i)+30) = zeros(1,31);
    end
end

% FOR PLOTTING
% Compare the new data with raw data
A = 1; B = 3e6;
t = (A:B)/3e1;
ch = 8;
% plot(t,raw(ch,A:B));
% hold on;
plot(t, detrended(ch,A:B)); 
hold on;
plot(t, blanked(ch,A:B))
xlabel('Time (ms)'); ylabel('Amplitude (uV)'); 
gridxy(dsc_times/3e1)


save(fullfile(save_path,'detrended.mat'),'detrended','-v7.3');
 %% 4. Blank 
 
load('sorted_TS_current.mat');
stim_times = sort(sorted_TS_current(:,1));

blanked = detrended;
for ch = 1:size(detrended,1)
    ch_data = detrended(ch,:);
    for i = 1:numel(stim_times)
       blanked(ch,stim_times(i):stim_times(i)+30) = zeros(1,31);
    end
    
    for j = 1:numel(dsc_times)
       blanked(ch,dsc_times(j)-3:dsc_times(j)+3) = zeros(1,7);
    end
end

%%
% clear new_data

save_path = uigetdir;
% save(fullfile(save_path,'detrended_data.mat'), 'copy', '-v7.3')
writemda32(blanked,fullfile(save_path,'preprocessed_data.mda'));