function plotStimPulses(raw,stimInfo,arrs,channelInfo)
    
    stimChannelInfo = channelInfo.singleStimCh;
    corrected_stim_idx = stimChannelInfo.singleRowIdx;
    corrected_rec_idx = channelInfo.allCh.allRowIdx;
    singleChannelNames = stimChannelInfo.singleChNames;
    stimPulseArr = arrs.singleStimPulseArr;
    stim_pulses = stimInfo(2,:);
    
    nCH = size(stimPulseArr,1); nCU = size(stimPulseArr,2);
    nStim = numel(stim_pulses);
    % 20 
    stim_ch_idx = 20; % index with respect to stim channels
    real_stim_ch = corrected_stim_idx(stim_ch_idx);
    rec_ch_idx = real_stim_ch; % index with respect to all channels 
    rec_ch_idx = 11
    currents = [0.5,1,2,5];
    T1 = 0; T2 = 800; 
    for i = 1:4
        figure;
        subplot(2,1,1)
        ts = cell2mat(stimPulseArr(stim_ch_idx,i));
        raw_trials = zeros(10,T1+T2+1);
        for j = 1:numel(ts)
            raw_trials(j,:) = raw(rec_ch_idx,ts(j)-T1:ts(j)+T2);
        end
        mean_sub = raw_trials-mean(raw_trials);
%         plot(mean_sub')
        plot(raw_trials')
        title(sprintf('10 reps of stim channel %s and current %.1f',singleChannelNames{stim_ch_idx},currents(i)))
        xlabel('Samples')
        ylabel('uV')
        subplot(2,1,1)
        dx = gradient(mean_sub,1);
%         plot(dx')
%         title(sprintf('Derivatives of each trial of stim channel %s and current %.1f',singleChannelNames{stim_ch_idx},currents(i)))
%         xlabel('Samples')
%         ylabel('uV')
        
%         subplot(2,1,2)
        y = movvar(dx',30);
        a = find(mean(y') < 80); a = a(2);
        
%         semilogy(y)
%         xline(a,'LineWidth',1.5)
%         hold on;
%         title('Moving variance of window size 30 (log scale) with proposed cutoff point using threshold')
%         subplot(2,1,1)
        blanked1 = [zeros(10,a-2) mean_sub(:,a-1:end)];
        blanked2 = mean_sub(:,a-1:end);
%         plot(blanked1')
%         plot(blanked1(trial,:)')
%         title(sprintf('Blanked trial of stim channel %s and current %.1f',singleChannelNames{stim_ch_idx},currents(i)))
        subplot(2,1,2)

        dtrend = [];
        for k = 1:10
            b = blanked2(k,:);
            db = nt_detrend(b',5,[],[],[],[],60)';
            z = 0;
            while ~((db(1) > -25))
                db(1) = [];
                z = z + 1;
            end
            dtrend(k,:) = [zeros(1,a-2) zeros(1,z) db];
        end
        plot(dtrend')
% %         plot(dtrend(trial,:)')
%         title('Detrended data')
    end
        
end
    

    