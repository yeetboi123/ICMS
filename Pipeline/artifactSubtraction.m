%% Add paths to directories
addpath('mdaio/')
addpath('neuroshare/')
addpath('basic_util/')
save_path = uigetdir('Save directory')
%% Load recorded stimulation data and stimulation times
[raw_data, ~, stim_data, ~]=read_Intan_RHS2000_file;
% Assume Intan data already loaded in parent script 
nCH = size(raw_data,1);
nSamples = size(raw_data,2);
STIM_PER_BLOCK = 1;
NUM_REPLICANTS = 11;
nCU = 4;
currents = [0.5 1 2 5];

[time_arr,current_arr] = getCathTimeNew(stim_data); time_arr_s = time_arr/3e4;
time_vector_unsorted = reshape(time_arr', 1, []); current_vector = reshape(current_arr', 1, []);
temp = repmat(new_order',1,nCU*NUM_REPLICANTS)'; channel_vector = temp(:)';

[B,current] = sort(time_vector_unsorted','ascend');
sorted_TS_current = [B current_vector(current)' channel_vector(current)'];
setGlobalx(sorted_TS_current)

time_vector = sort(time_arr(:));
time_vector_s = time_vector/3e4;
time_vector_ms = time_vector/3e1;
%% Plot stimulation pulses and label with stimulating channel name 
subplot(2,1,1)
plotSamples = 360000;
plot((1:plotSamples)/3e4,raw_data(:,1:plotSamples)');
% I need to get channel names somehow... 
new_order

for ch = 1:nCH
    xline(time_arr_s(ch,1),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',new_order(ch),current_arr(ch,1))); 
    hold on;
    xline(time_arr_s(ch,2),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',new_order(ch),current_arr(ch,2))); 
    hold on;
    xline(time_arr_s(ch,3),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',new_order(ch),current_arr(ch,3))); 
    hold on;
    xline(time_arr_s(ch,4),'Color','r','Linestyle','--','Linewidth',2,...
        'label',sprintf('Ch.%d @ %.1f uA',new_order(ch),current_arr(ch,4))); 
    hold on;
end
title('Amplifier data for ICMS 16 for 4 replicants with stimulation channel and current')
xlim([0 plotSamples/3e4])

subplot(2,1,2)
plot((1:plotSamples)/3e4,stim_data(:,1:plotSamples)');
xlabel('Time (s)')
ylabel('Amplitude (uV)')
title('Randomized stimulation pulse amplitudes for ICMS 16 for 4 replicants')
xlim([0 plotSamples/3e4])

%% Plot all channels on shank 
% Uncomment to allow plotting
total_samples = 120000;
T = double(1:120000)/30;
% plot(T, raw_data(:,1:total_samples)');
% hold on;
% 
% gridxy(stim_times_ms,'Color','r','Linestyle','--','Linewidth',2);
% xlabel('Time(ms)'); ylabel('Amplitude(uV)'); 
% title(sprintf('Plotting shank %d channels for 4 second window...',shank_index))

% Prompt user for blanking segment length
user_satisfied_flag = 'N';
while true
    
    if user_satisfied_flag ~= 'N'
        break
    else
        
        blank = input('How many ms to be blanked after stimulation? ');
        close all;
        new_data = raw_data(:,1:total_samples);
        for i = 1:100
            new_data(:,time_vector(i):time_vector(i)+double(blank)*30) = zeros(nCH,double(blank)*30+1);
        end
        plot(T, new_data(:,1:total_samples)');
        hold on;
        gridxy(time_vector_ms,'Color','r','Linestyle','--','Linewidth',2);
        
 
        xlabel('Time(ms)'); ylabel('Amplitude(uV)'); 
        title('Plotting 4 second window...')
        user_satisfied_flag = input('Is this OK? (Y/N): ','s');
        
    end
end


%% Separate windows based on stim channel and current 

% Grab segments around stimulation pulses
% Just midpoint stimulation pulse index
N = floor((length(time_vector)/STIM_PER_BLOCK)/2);
% Minimum length of window between stim pulse and subsequent stim pulse
% to assure shorter than all windows
stim_window = round(min(diff(time_vector)));
stim_window_last = nSamples - time_vector(end)-1;

% 4 dimensions: stim_combination (stim_ch + current combination), recording
% channel, samples, replicants
getCurrentInd = @(x) find(abs(sorted_TS_current(:,2) - x) < 0.001);
% Useful anon functions
getChInd = @(x) find(sorted_TS_current(:,3) == x);
window_array4 = zeros(nCH*nCU,nCH,stim_window+1,NUM_REPLICANTS);

ij_counter = 0;
% ij_counter specifies order of windows (e.g. first dimension of 1 is stim ch. 1 at 0.5 uA, 2 is stim ch. 1
% at 1 uA, 3 is stim ch. 1 at 2 uA, ...)
for i = 1:nCH
    for j = 1:nCU
        ij_counter = ij_counter + 1;
        stim_ch = new_order(i); stim_current = currents(j);
        stim_combo_ind = intersect(getCurrentInd(stim_current),getChInd(stim_ch));
        for k = 1:length(stim_combo_ind)
            stim_combo_idx = sorted_TS_current(stim_combo_ind(k),1);
            if stim_combo_idx == time_vector(end)
                window_array4(ij_counter,:,:,k) = [raw_data(:,stim_combo_idx:stim_combo_idx + stim_window_last) zeros(nCH,stim_window - stim_window_last)];
            else
                window_array4(ij_counter,:,:,k) = raw_data(:,stim_combo_idx:stim_combo_idx + stim_window);           
            end

            window_array4(ij_counter,:,1:double(blank)*30,k) = zeros(nCH,double(blank)*30);

        end
    end
end
%% Plot windows
figure;
plot_length = 10000;
rec_ch = 17;
stim_ch = 2;
for stim_combo_ind = [1 2 3 4] + ((stim_ch-1)*nCU)
    current_idx = stim_combo_ind - (((stim_ch-1)*nCU));
    current = currents(current_idx);
    subplot(4,1,current_idx)
    plot((1:plot_length)/3e1,squeeze(window_array4(stim_combo_ind,rec_ch,1:plot_length,:)))
    title(sprintf('Stim current: %.1f',current))
end
sgtitle(sprintf('%d replicants for recording ch. %d for stim ch. %d',NUM_REPLICANTS,rec_ch,stim_ch))
xlabel('Time (ms)')
ylabel('Amplitude (uV)')


%% Get median of windows 

disp('Calculating medianSignal...')
close all;

median_array = zeros(nCH*nCU,nCH,stim_window+1);

for stim_combo = 1:nCH*nCU
    for rec_ch = 1:nCH
        arr = squeeze(window_array4(stim_combo,rec_ch,:,:));
        median_array(stim_combo,rec_ch,:) = mean(arr,2);
    end
end

%% Plot median signals alongside raw data windows
figure;
plot_length = 10000;
rec_ch = 2;
stim_ch = 7;
for stim_combo_ind = [1 2 3 4] + ((stim_ch-1)*nCU)
    current_idx = stim_combo_ind - (((stim_ch-1)*nCU));
    current = currents(current_idx);
    subplot(4,1,current_idx)
    plot((1:plot_length)/3e1,squeeze(window_array4(stim_combo_ind,rec_ch,1:plot_length,:)))
    title(sprintf('Stim current: %.1f',current))
    hold on;
    plot((1:plot_length)/3e1,squeeze(median_array(stim_combo_ind,rec_ch,1:plot_length)),'r--','LineWidth',3)
end
sgtitle(sprintf('%d replicants for recording ch. %d for stim ch. %d',NUM_REPLICANTS,rec_ch,stim_ch))
xlabel('Time (ms)')
ylabel('Amplitude (uV)')

%% Plot just the median signals
figure;
plot(squeeze(median_array(1,1:19,:))')
title('Median')
%% PCA analysis 

% Window array is 3-dimensional (x - stim indices, y - time, z - channel)
X = [];
for i = 1:76
    for j = 1:19
        x = squeeze(median_array(i,j,50:300))';
        X = [X;x];
    end
end
labels = repmat(1:4,76*19/4,1)'; labels = labels(:);
[coeff,score,latent,~,explained] = pca(X);

scatter3(score(:,1),score(:,2),score(:,3),19,labels,'filled'); colormap(parula)
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
title('5 Millisecond Post-Stim Windows')

%% TSNE
X = [];
for i = 1:19
    x = window_array4(:,50:150,i);
    X = [X;x];
end
labels = repmat(1:19,2280,1); labels = labels(:);

Y3 = tsne(X,'Algorithm','barneshut','NumPCAComponents',20,'NumDimensions',3,'Perplexity',20);
figure
% gscatter(Y3(:,1),Y3(:,2),labels);

scatter3(Y3(:,1),Y3(:,2),Y3(:,3),19,labels,'filled');colormap(hot)
legend

%% Plot K_plot windows with mean of K_windows windows subtracted 

temp1 = zeros(nCH*nCU,nCH,size(window_array4,3),size(window_array4,4));
for stim_combo = 1:nCH*nCU
    parfor rec_ch = 1:nCH
        temp1(stim_combo,rec_ch,:,:) = window_array4(stim_combo,rec_ch,:,:) - median_array(stim_combo,rec_ch,:);
    end
end


figure;
plot_length = 21840;
rec_ch = 9;
stim_ch = 9;
for stim_combo_ind = [1 2 3 4] + ((stim_ch-1)*nCU)
    current_idx = stim_combo_ind - (((stim_ch-1)*nCU));
    current = currents(current_idx);
    
    subplot(4,2,current_idx)
    plot((1:plot_length)/3e1,squeeze(window_array4(stim_combo_ind,rec_ch,1:plot_length,:)))
    title(sprintf('Stim current: %.1f',current))
    hold on;
    plot((1:plot_length)/3e1,squeeze(median_array(stim_combo_ind,rec_ch,1:plot_length)),'r--','LineWidth',3)
end

for stim_combo_ind = [1 2 3 4] + ((stim_ch-1)*nCU)
    current_idx = stim_combo_ind - (((stim_ch-1)*nCU));
    current = currents(current_idx);
    
    subplot(4,2,current_idx + 4)
    plot((1:plot_length)/3e1,squeeze(temp1(stim_combo_ind,rec_ch,1:plot_length,:)))
    title(sprintf('Median subtracted at stim current: %.1f',current))
end

sgtitle(sprintf('%d median subtracted replicants for recording ch. %d for stim ch. %d',NUM_REPLICANTS,rec_ch,stim_ch))
xlabel('Time (ms)')
ylabel('Amplitude (uV)')

%% Process data for all windows

% Iterate through every stimulation pulse and replace window around it with
% new processed window
disp('Creating new data array...');

new_data = raw_data; 
for channel = 1:nCH
    for stim_ind = 1:numel(time_vector)
        [stim_combo_ind, replicant] = reverseMap(time_vector(stim_ind),new_order,nCU);
        new_data(channel,time_vector(stim_ind):time_vector(stim_ind)+stim_window) = ...
            temp1(stim_combo_ind,channel,:,replicant);
    end
%     new_data(channel,start_burst(end):end) = zeros(1,length(new_data)-time_vector(end)+1);
end

for channel = 1:10
    h = figure;
    plot((1:length(new_data))/3e4,new_data(channel,:)'); hold on; ...
        plot((1:length(new_data))/3e4,raw_data(channel,:)')
    legend('Artifact subtraction','Original')
    xlabel('Time(s)')
    ylabel('Amplitude(uV)')
    title(sprintf('Comparison between original and artifact-subtracted Ch. %d data', channel))
end
uiwait(h);
%% Convert into .mda file
disp('Saving new data as mda file...');
% Select specific channels
writemda32(new_data, fullfile(save_path,'preprocessed_data.mda'));
save(fullfile(save_path,'sorted_TS_current.mat'),'sorted_TS_current')
%% Done!
disp('Done!');


%% Helper functions

function [stim_combo_ind, replicant] = reverseMap(time,new_order,nCU)
%     global sorted_TS_current
    sorted_TS_current = getGlobalx;
    idx = find(sorted_TS_current(:,1) == time);
    current = sorted_TS_current(idx, 2);
    stim_ch = sorted_TS_current(idx, 3);
    
    switch round(current,1)
        case 0.5
            offset = 1;
        case 1
            offset = 2;
        case 2 
            offset = 3;
        case 5
            offset = 4;
    end

    rel_ch = find(new_order == stim_ch);
    stim_combo_ind = offset + (rel_ch-1)*nCU;
    
    idx2 = find((abs(sorted_TS_current(:,2) - current) < 0.001) & (sorted_TS_current(:,3) == stim_ch));
    replicant = find(idx2 == idx);
end

function setGlobalx(val)
    global x
    x = val;
end

function r = getGlobalx
    global x
    r = x;
end