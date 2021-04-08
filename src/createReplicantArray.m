function replicant_array = createReplicantArray(raw,stimInfo,channelInfo,stimPulseArrs,nSamples,debug_flag)
    % create input array for artifact removal 
    
    stimChannelInfo = channelInfo.singleStimCh;
    corrected_stim_idx = stimChannelInfo.singleRowIdx;
    corrected_rec_idx = channelInfo.allCh.allRowIdx;
    singleChannelNames = stimChannelInfo.singleChNames;
    
    stimPulseArr = stimPulseArrs.singleStimPulseArr; 
    stim_pulses = stimInfo(2,:);
    
    currents = [0.5,1,2,5];
    nSCH = size(stimPulseArr,1); nCU = numel(currents);
    nRCH = size(raw,1); nRep = numel(stimPulseArr{1,1});
    
    replicant_array = zeros(nSCH,nRCH,nCU,nRep,nSamples);
    for i = 1:nSCH
        fprintf('%d stim channel/%d total stim channels...\n',i,nSCH)
        for j = 1:nRCH
            for k = 1:nCU
                stim_times = stimPulseArr{i,k};
                for l = 1:nRep
                    replicant_array(i,j,k,l,:) = raw(j,stim_times(l):stim_times(l)+nSamples-1);        
                end
            end
        end
    end
    
    if debug_flag
        debug(20,10);
    end
    
    function debug(stim_ch_idx,rec_ch_idx)
        if nargin == 1
            rec_ch_idx = corrected_stim_idx(stim_ch_idx);
        end
        close all
%         T1 = 0; T2 = 800; 
        rec_ch_idx = corrected_rec_idx(rec_ch_idx);
        for kk = 1:4
            figure;         
            ts = cell2mat(stimPulseArr(stim_ch_idx,kk));
            raw_trials = squeeze(replicant_array(stim_ch_idx,rec_ch_idx,kk,:,:));
            plot(raw_trials')
            title(sprintf('10 reps of stim channel %s and current %.1f',singleChannelNames{stim_ch_idx},currents(kk)))
            xlabel('Samples')
            ylabel('uV')
        end
    end

end

% input_data = createInputData(raw,stimInfo,channelInfo,stimPulseArr,0);

