% Artifact mitigation script

% 1. Raw data 
% 2. Mean subtraction to get rid of movement artifacts
% 3. nt_detrend with degree 3 polynomial at SD of 6, 3 iterations, 60 sample window 
% 4. Blank 3 ms before and 3.5 ms after and then interpolate using pchip
% 5. Put into Mountainsort 

%% 1. Import raw data 

 [raw, ~, ~, ~]=read_Intan_RHS2000_file;
 clear amplifier_channels amp_settle_data charge_recovery_data compliance_limit_data stim_data t 
 
 %%  2. Channel mean subtraction 
 
 mean_sg = mean(raw);
 raw = raw - mean_sg; clear mean_sg
 
 %% 3. Detrend (takes a while)
 
NUM_CORES = 12;
L = numel(raw(1,:));
segLength = floor(L/NUM_CORES);
new_data = zeros(size(raw,1),L);

for CH_IND = 1:19
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
        segArr2{seg} = nt_detrend(segArr1{seg}',3,[],[],6,3,60);
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

 %% 4. Blank and interpolate 
 
 load('sorted_TS_current.mat');
 timestamps = sort(sorted_TS_current(:,1));
 

copy = new_data;
preblank = 3; postblank = 3.5;

for ch = 1;size(new_data,1);
    for stim_idx = 1:length(timestamps) - 1
        x = 1:16*3e1+1;
        prewindow = copy(ch,timestamps(stim_idx)-8*3e1: timestamps(stim_idx)+8*3e1);
        % blank before and after stim pulse
        copy(ch,timestamps(stim_idx)-preblank*3e1: timestamps(stim_idx)+postblank*3e1) = zeros(1,(preblank+postblank) *3e1 + 1);
        % interpolate samples in the window 
        window = copy(ch,timestamps(stim_idx)-8*3e1: timestamps(stim_idx)+8*3e1);
        idx = window~=0;

        yn = interp1(x(idx),window(idx),x,'pchip');
        copy(ch,timestamps(stim_idx)-8*3e1: timestamps(stim_idx)+8*3e1) = yn;
    end
end
 
clear new_data

save_path = uigetdir;
save(fullfile(save_path,'detrended_data.mat'), 'copy')