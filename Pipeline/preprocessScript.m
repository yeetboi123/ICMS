% Artifact mitigation script

%% 1. Import raw data 
addpath mdaio

[raw, ~, ~, ~]=read_Intan_RHS2000_file;
clear amplifier_channels amp_settle_data charge_recovery_data compliance_limit_data stim_data t 
 
 %%  2. Channel mean subtraction (for noisy data)
mean_sg = mean(raw);
med_sg = median(raw);
raw = raw - mean_sg; clear mean_sg
t = (12e4:24e4)/3e1;
subplot 211
plot(t,raw(8,12e4:24e4)); xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Raw data')
subplot 212
plot(t,raw(8,12e4:24e4)- mean_sg(12e4:24e4)); xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Mean-subtracted data')
sgtitle('ICMS15: Raw data vs. mean-subtracted data')

subplot 211
plot(raw(:,1:30e4)')
subplot 212
plot((raw(:,1:30e4) - med_sg(1:30e4))')

%% 2a. Remove strange discontinuities seen in Intan recordings 
data = raw(1,1:3600000);
for stim_idx = 1:8
    start_idx = stim_times(stim_idx) + 100;
    end_idx = stim_times(stim_idx + 1) - 100;
    y = abs(diff(data(start_idx:end_idx)));
    [~,ind] = max(y);
    corrected_ind = ind + start_idx;
    data(corrected_ind-50:corrected_ind+50) = zeros(1,101);
end
plot(raw(1,1:3600000));hold on; plot(data)
plot(data(start_idx:end_idx)); hold on
plot(data(1:360000))

plot(diff(data)); hold on; plot(data)

%% 3. Detrend (takes a while)
addpath NoiseTools\
 
NUM_CORES = 12;
L = numel(raw(1,:));
segLength = floor(L/NUM_CORES);
new_data = zeros(size(raw,1),L);

for CH_IND = 1:size(raw,1)
    CH_IND 
    tic
    segArr1 = cell(12,1);
    segArr2 = cell(12,1);
    ch = raw(CH_IND,:);
    
    for seg = 1:NUM_CORES
        startInd = segLength * (seg - 1) + 1;
        endInd = startInd + segLength - 1;

        if seg == NUM_CORES
            segArr1{seg} = ch(startInd:end);
        else
            segArr1{seg} = ch(startInd:endInd);
        end
    end

    parfor seg = 1:NUM_CORES
        seg
        tic
        segArr2{seg} = nt_detrend(segArr1{seg}',5,[],[],6,3,60);
        toc
    end

    reconstruct = [];
    for seg = 1:NUM_CORES
        reconstruct = [reconstruct segArr2{seg}'];
    end
    
    new_data(CH_IND,:) = reconstruct; clear reconstruct
    toc
end
 
clear raw;
% 
% subplot 211
% t = (1:6e4)/3e1;
% plot(t,raw(1,1:6e4));xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Raw data')
% xlim([800 1100]);
% ylim([-1500 3000])
% subplot 212
% detrended_data = nt_detrend(ch1',5,[],[],6,3,60);
% plot(t,detrended_data);xlabel('Time (ms)'); ylabel('Amplitude (uV)'); title('Detrended data')
% xlim([800 1100]);
% ylim([-1000 1000])
% sgtitle('ICMS15: Raw data vs. detrended data')
 %% 4. Blank and interpolate 
 
load('sorted_TS_current.mat');
timestamps = sort(sorted_TS_current(:,1));

window = 120;
stim_buffer = 30;
copy = new_data;

% for ch = 1:size(copy,1)
    
good_spikes = []
for ch = 3
    
%     for stim_idx = 1:length(timestamps) - 1
    for stim_idx = 1:30
        ts = timestamps(stim_idx);
        % detrend segment before and after
        prestimSeg = copy(ch,ts-window:ts-stim_buffer/2);
        poststimSeg = copy(ch,ts+stim_buffer:ts+window);
        
        [prestimSegDT,~] = nt_detrend(prestimSeg',5,[],[],[],[],[]);
        [poststimSegDT,~] = nt_detrend(poststimSeg',5,[],[],[],[],[]);

        copy(ch,ts-window:ts-stim_buffer/2) = prestimSegDT;
        copy(ch,ts+stim_buffer:ts+window) = poststimSegDT;

        [istartPre1,istopPre1,distPre1] = findsignal(prestimSegDT,spike_template,'TimeAlignment','dtw','Metric','euclidean','MaxNumSegments',2);
        [istartPre2,istopPre2,distPre2] = findsignal(prestimSegDT,spike_template2,'TimeAlignment','dtw','Metric','euclidean','MaxNumSegments',2);
        [istartPost1,istopPost1,distPost1] = findsignal(poststimSegDT,spike_template,'TimeAlignment','dtw','Metric','euclidean','MaxNumSegments',2);
        [istartPost2,istopPost2,distPost2] = findsignal(poststimSegDT,spike_template2,'TimeAlignment','dtw','Metric','euclidean','MaxNumSegments',2);
        
        distPre = [distPre1 distPre2]; istopPre = [istopPre1 istopPre2]; 
        distPost = [distPost1 distPost2]; istartPost = [istartPost1 istartPost2];
        
        [idxPre,good_distPre] = find(distPre < 450);
        [idxPost,good_distPost] = find(distPost < 450);
        if size(good_distPre,1) > 0
            seg_start = window-max(istopPre(idxPre));
            good_spikes = [good_spikes ts - seg_start];
            copy(ch,ts - seg_start:ts) = zeros(1,seg_start+1);
        else
            copy(ch,ts-window:ts) = zeros(1,window+1);
        end
        
        if size(good_distPost,1) > 0
            seg_end = min(istartPost(idxPost)) + stim_buffer;
            good_spikes = [good_spikes ts + seg_end];
            copy(ch,ts:ts + seg_end) = zeros(1,seg_end+1);
        else
            copy(ch,ts:ts+window) = zeros(1,window+1);
        end

    end
end
plot(new_data(3,1:1e6)'); hold on
% figure
plot(copy(3,1:1e6),'Linewidth',1);
gridxy(timestamps)
gridxy(timestamps + stim_buffer)
gridxy(timestamps - stim_buffer/2)
gridxy(good_spikes,'Color','r')




x1 = 60e4; x2 =120e4; ch = 1:15;
t = (x1:x2)/30;
subplot 211
plot(t,new_data(ch,x1:x2)');
subplot 212
plot(t,copy(ch,x1:x2)');
% legend({'Detrended','Interpolated'})
xlabel('Time (ms)')
ylabel('Amplitude (uV)');
title('ICMS15: Detrended vs. Blanked and Interpolated')
gridxy(timestamps/30)

figure;

%%
% clear new_data

save_path = uigetdir;
% save(fullfile(save_path,'detrended_data.mat'), 'copy', '-v7.3')
writemda32(copy,fullfile(save_path,'preprocessed_data.mda'));