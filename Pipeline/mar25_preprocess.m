addpath basic_util/ 
addpath ICMS15/Mar24_ICMS15/Pipeline
addpath mdaio/
[raw,stim_data,new_order,stim_new_order,t,save_path] = createMappingNEW;
clear amp_settle_data charge_recovery_data compliance_limit_data 

subplot 211
plot(raw(1:20,1:1e6)')
gridxy(find(t==0));
subplot 212
stim_chs=find(any(stim_data,2) == 1);
plot(stim_data(stim_chs,1:1e6)')
gridxy(find(t==0));
%% Get stim times
[time_arr,current_arr] = getCathTimeNew(stim_data);
stim_times = sort(reshape(time_arr', 1, []));

nCH = size(raw,1);
nSamples = size(raw,2);
STIM_PER_BLOCK = 1;
NUM_REPLICANTS = 10;
nCU = 4;
currents = [0.5 1 2 5];

time_arr_s = time_arr/3e4;
time_vector_unsorted = reshape(time_arr', 1, []); current_vector = reshape(current_arr', 1, []);
temp = repmat(stim_new_order',1,nCU*NUM_REPLICANTS)'; channel_vector = temp(:)';

[B,current] = sort(time_vector_unsorted','ascend');
sorted_TS_current = [B current_vector(current)' channel_vector(current)'];
save(fullfile(save_path,'sorted_TS_current.mat'),'sorted_TS_current')

%%
subplot(2,1,1)
plotSamples = 360000;
plot((1:plotSamples)/3e4,raw(:,1:plotSamples)');
% I need to get channel names somehow... 

for ch = 1:numel(stim_new_order)
    xline(time_arr_s(ch,1),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',stim_new_order(ch),current_arr(ch,1))); 
    hold on;
    xline(time_arr_s(ch,2),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',stim_new_order(ch),current_arr(ch,2))); 
    hold on;
    xline(time_arr_s(ch,3),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',stim_new_order(ch),current_arr(ch,3))); 
    hold on;
    xline(time_arr_s(ch,4),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',stim_new_order(ch),current_arr(ch,4))); 
    hold on;
end
title('Amplifier data for ICMS 16 for 4 replicants with stimulation channel and current')
xlim([0 plotSamples/3e4])

subplot(2,1,2)
plot((1:plotSamples)/3e4,stim_data(:,1:plotSamples)');

% for k = 1:22
%     xline(stim_times(k)/3e4,'Color','r','Linestyle','--','Linewidth',2,...
%         'label',sprintf('Ch.%d',new_order(stim_order(k)))); 
%     hold on;
% end

chs = sorted_TS_current(:,3);
for k = 1:22
    xline(stim_times(k)/3e4,'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d',chs(k))); 
    hold on;
end


xlabel('Time (s)')
ylabel('Amplitude (uV)')
title('Randomized stimulation pulse amplitudes for ICMS 16 for 4 replicants')
xlim([0 plotSamples/3e4])
%% Test median/mean subtraction
mean_sig = mean(raw);
meanSub = raw - mean_sig;
plot(meanSub(:,1:2e5)')
figure
plot(raw(:,1:2e5)')

plot(mean_sig(1:2e5))
%%
addpath NoiseTools
NUM_CORES = 12;
L = length(stim_times);
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
% for CH_IND = 1
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
            data_arr{c} = sample(stim_times(t1):end); % last segment stretches to the end
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
        norm_stim_times = stim_times(sIDX) - stim_times(sIDX(1)) + 1; % normalize
%         if seg == NUM_CORES
%             sIDX(end) = [];
%         end
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
                y1 = nt_detrend(chunk1',5,[],[],[],[],60);
                ys = [ys zeros(1,30) y1'];
            else
                chunk1 = segment(norm_stim_times(i)+30:norm_dsc_times(i)-1); 
                chunk2 = segment(norm_dsc_times(i):norm_stim_times(i+1)-1);
                y1 = nt_detrend(chunk1',5,[],[],[],[],60);
                y2 = nt_detrend(chunk2',5,[],[],[],[],60);
                y2(1:15) = ones(1,15);
                ys = [ys zeros(1,30) y1' y2'];
            end
            
        end
        yeet = [yeet ys];
    end
%     end_piece = sample(stim_times(end):end);
%     end_piece_dt = nt_detrend(end_piece',5,[],[],[],[],60);
    detrended(CH_IND,:) = [zeros(1,stim_times(1)-1) yeet];
    toc
end

x2 = 3e5;
plot(detrended(1,1:x2))
hold on;
plot(raw(1,1:x2))
gridxy(stim_times)


save(fullfile(save_path,'detrended_no_mean_sub.mat'),'detrended_no_mean_sub','-v7.3');
%% Test SVM 
load('SVMmodel.mat')
stim_buffer = 40;
threshold = 0.5;
ch = 1;
for i = 1:100
    
seg_raw = raw(ch,stim_times(i)+stim_buffer:stim_times(i)+stim_buffer+30*5);    
seg = detrended(ch,stim_times(i)+stim_buffer:stim_times(i)+stim_buffer+30*5);
y = nt_detrend(seg',5,[],[],[],[],60);
yr = nt_detrend(seg_raw',5,[],[],[],[],60);
vec = y;
vec_og = seg;
vec_r = yr;
X_test = createRollingWindow(vec,36);
X_verify = createRollingWindow(vec_og,36);
Xr_test = createRollingWindow(vec_r,36);

trainedSVM = trainedModel.ClassificationSVM;

[labels,scores] = predict(trainedSVM, X_test);
[labelsR,scoresR] = predict(trainedSVM, Xr_test);
% ad hoc conditions
labelCorrect = labels == 1;
passesThreshold = scores(:,2) > threshold;
goodAmplitude = min(X_test') > -480;
goodAmplitudeOG = range(X_verify') < 300;
good_idx = find((labelCorrect & passesThreshold & goodAmplitude' & goodAmplitudeOG') == 1);
 
all_segments = [good_idx good_idx + 36];

good_scores = scores(good_idx);

earliest_segment = [min(good_idx) min(good_idx)+36];
scores1 = scores(:,2);
earliest_seg_score = scores1(min(good_idx));


set(gcf,'Visible','on')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

plot(vec_og)
hold on;
plot(vec); 
title(sprintf('Max score: %.1f',earliest_seg_score))
gridxy(earliest_segment,'Color','r')
gridxy(all_segments); 
uiwait(gcf);
end
%% 
load('SVMmodel.mat')
nCH = size(detrended,1);
SVM_corrected = detrended;
threshold = 0.5;
stim_buffer = 40;

ch_with_spike = zeros(numel(stim_times),1);
scores_arr = zeros(numel(stim_times),1);

for i = 1:numel(stim_times)
% for i = 1:5
  stim_blank_idx = zeros(nCH,numel(stim_times));
  scores2 = zeros(nCH,numel(stim_times));
  stim_time = stim_times(i);
  
  segs = detrended(:,stim_time+stim_buffer:stim_time+stim_buffer+30*5);
%   segs_raw = raw(:,stim_time+stim_buffer:stim_time+stim_buffer+30*5);
  
  for j = 1:nCH
    seg = segs(j,:);
%     raw_seg = segs_raw(j,:);
    y = nt_detrend(seg',6,[],[],[],[],60);
    windows = createRollingWindow(y,36);
    X_verify = createRollingWindow(seg,36);
    trainedSVM = trainedModel.ClassificationSVM;
    [labels,scores] = predict(trainedSVM, windows);

    labelCorrect = labels == 1;
    passesThreshold = scores(:,2) > threshold;
    goodAmplitude = min(windows') > -480;
%     goodAmplitude2 = max(windows') < 200;
    goodRangeOG = range(X_verify') < 300;
    scores = scores(:,2);
    good_idx = find((labelCorrect & passesThreshold & goodAmplitude' & goodRangeOG') == 1);
    good_scores = scores(good_idx); 

    [stop_blank,idx2] = min(good_idx);

    if ~isempty(stop_blank) 
      stim_blank_idx(j,i) = stop_blank;
      scores2(j,i) = good_scores(idx2);
    else
      stim_blank_idx(j,i) = 121; 
      scores2(j,i) = 0;
    end
  end
  
  [stop_blank_min,idx] = min(stim_blank_idx(:,i));
  ch_with_spike(i) = idx;
  scores_arr(i) = scores2(idx,i);
  SVM_corrected(:,stim_times(i):stim_times(i)+stim_buffer+stop_blank_min-5) =...
      zeros(nCH,stim_buffer+stop_blank_min-5+1);

end 
%%
% ch = 16:32
% x1 = 1
% x2 = 1e5
% t = (x1:x2)/3e1;
% plot(t,detrended(ch,x1:x2)+300,'LineWidth',2);
% hold on;
% plot(t,SVM_corrected(ch,x1:x2));
% text(stim_times/3e1,zeros(numel(stim_times),1),string(scores_arr));
% text(stim_times/3e1,zeros(numel(stim_times),1)+100,string(ch_with_spike));
% % hold on;
% % plot(t,raw(ch,1:x2)-200);
% gridxy(stim_times/3e1)
% xlabel('Time (ms)')
% ylabel('Amplitude (uV)')
% title('Selective blanking w/ SVM example')


for k = 1:100
    ch = ch_with_spike(k);
    x1 = stim_times(k)-500; x2 = stim_times(k)+500;
    t = (x1:x2)/3e1;
    plot(t,detrended(ch:ch+3,x1:x2)'+200,'LineWidth',1);
    hold on;
    plot(t,SVM_corrected(ch:ch+3,x1:x2)');
    text(stim_times(k)/3e1,50,string(scores_arr(k)));
    xline(stim_times(k)/3e1)
    uiwait(gcf);
end

writemda32(SVM_corrected,fullfile(save_path,'preprocessed_data_single.mda'));
%% Preprocess paired stim data
[raw, amplifier_channels, stim_data, t]=read_Intan_RHS2000_file;
clear amp_settle_data charge_recovery_data compliance_limit_data
[time_arr,current_arr] = getCathTimePaired(stim_data);
stim_times = unique(cell2mat(time_arr));

subplot 211
plot(raw(1:20,1:3e5)')
gridxy(find(t==0));
gridxy(stim_times,'color','r')
subplot 212
stim_chs=find(any(stim_data,2) == 1);
for i = 1:numel(stim_chs)
    plot_ch = stim_chs(i);
    plot(stim_data(plot_ch,1:3e5)'+i*6,'LineWidth',3);
    hold on
end
gridxy(find(t==0));
gridxy(stim_times,'color','r')


%% Test median/mean subtraction
mean_sig = mean(raw);
meanSub = raw - mean_sig;
plot(meanSub(:,1:2e5)')
figure
plot(raw(:,1:2e5)')
%%
addpath NoiseTools\
NUM_CORES = 12;
L = length(stim_times);
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
% for CH_IND = 1
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
            data_arr{c} = sample(stim_times(t1):end); % last segment stretches to the end
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
        norm_stim_times = stim_times(sIDX) - stim_times(sIDX(1)) + 1; % normalize
%         if seg == NUM_CORES
%             sIDX(end) = [];
%         end
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
                y1 = nt_detrend(chunk1',5,[],[],[],[],60);
                ys = [ys zeros(1,30) y1'];
            else
                chunk1 = segment(norm_stim_times(i)+30:norm_dsc_times(i)-1); 
                chunk2 = segment(norm_dsc_times(i):norm_stim_times(i+1)-1);
                y1 = nt_detrend(chunk1',5,[],[],[],[],60);
                y2 = nt_detrend(chunk2',5,[],[],[],[],60);
                y2(1:15) = ones(1,15);
                ys = [ys zeros(1,30) y1' y2'];
            end
            
        end
        yeet = [yeet ys];
    end
%     end_piece = sample(stim_times(end):end);
%     end_piece_dt = nt_detrend(end_piece',5,[],[],[],[],60);
    detrended(CH_IND,:) = [zeros(1,stim_times(1)-1) yeet];
    toc
end

x2 = 3e5;
plot(detrended(1,1:x2))
hold on;
plot(raw(1,1:x2))
gridxy(stim_times)


save(fullfile(save_path,'detrended_paired.mat'),'detrended','-v7.3');

%% 

nCH = size(detrended,1);
SVM_corrected = detrended;
threshold = 0.8;
stim_buffer = 40;

ch_with_spike = zeros(numel(stim_times),1);
scores_arr = zeros(numel(stim_times),1);

for i = 1:numel(stim_times)
% for i = 1:5
  stim_blank_idx = zeros(nCH,numel(stim_times));
  scores2 = zeros(nCH,numel(stim_times));
  stim_time = stim_times(i);
  
  segs = detrended(:,stim_time+stim_buffer:stim_time+stim_buffer+30*5);
%   segs_raw = raw(:,stim_time+stim_buffer:stim_time+stim_buffer+30*5);
  
  for j = 1:nCH
    seg = segs(j,:);
%     raw_seg = segs_raw(j,:);
    y = nt_detrend(seg',6,[],[],[],[],60);
    windows = createRollingWindow(y,36);
    trainedSVM = trainedModel.ClassificationSVM;
    [labels,scores] = predict(trainedSVM, windows);

    labelCorrect = labels == 1;
    passesThreshold = scores(:,2) > threshold;
    goodAmplitude = min(abs(windows)') > -480;
    goodAmplitude2 = max(windows') < 200;
    scores = scores(:,2);
    good_idx = find((labelCorrect & passesThreshold & goodAmplitude' & goodAmplitude2') == 1);
    good_scores = scores(good_idx); 

    [stop_blank,idx2] = min(good_idx);

    if ~isempty(stop_blank) 
      stim_blank_idx(j,i) = stop_blank;
      scores2(j,i) = good_scores(idx2);
    else
      stim_blank_idx(j,i) = 121; 
      scores2(j,i) = 0;
    end
  end
  
  [stop_blank_min,idx] = min(stim_blank_idx(:,i));
  ch_with_spike(i) = idx;
  scores_arr(i) = scores2(idx,i);
  SVM_corrected(:,stim_times(i):stim_times(i)+stim_buffer+stop_blank_min-5) =...
      zeros(nCH,stim_buffer+stop_blank_min-5+1);

end 
writemda32(SVM_corrected,fullfile(save_path,'preprocessed_data_paired_only.mda'));
%%
paired_stim_ch_arr = [];
current_arr = [];
time_arr = [];
[time_arr,current_arr] = getCathTimePaired(stim_data)

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


for i = 1:numel(stim_times) 
    a = []; b = []; c = [];
    for ch = 1:numel(cathodes)
        if stim_data(cathodes(ch),stim_times(i)) < 0
            a = [a cathodes(ch)];
            b = [b stim_data(cathodes(ch),stim_times(i))];
            c = [c stim_times(i)];
        end
    end   
    paired_stim_ch_arr(i,:) = a;
    current_arr(i,1) = b(1); 
    time_arr(i,1) = c(1);
end

coords = readmatrix('channel_maps/128channel4shank_location_data.csv');
coords_of_interest = coords(new_order,:);
num2cell(coords_of_interest(stim_chs,:)

% 1st col = time, 2nd col = current, 3rd,4th = ch 1-idx, 5th, 6th = ch
% 128-order idx
sorted_TS_current_paired = [stim_times' current_arr paired_stim_ch_arr new_order(paired_stim_ch_arr)];
save(fullfile(save_path,'sorted_TS_current_paired.mat'),'sorted_TS_current_paired')

%% Combine

writemda32(preprocessed_data_combined,fullfile(save_path,'preprocessed_data_combined.mda'));
%%
function output = createRollingWindow(vector, n)
% CREATEROLLINGWINDOW returns successive overlapping windows onto a vector
%   OUTPUT = CREATEROLLINGWINDOW(VECTOR, N) takes a numerical vector VECTOR
%   and a positive integer scalar N. The result OUTPUT is an MxN matrix,
%   where M = length(VECTOR)-N+1. The I'th row of OUTPUT contains
%   VECTOR(I:I+N-1).
l = length(vector);
m = l - n + 1;
output = vector(hankel(1:m, m:l));
end