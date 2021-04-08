function [detrended_array,output_data,nSamples] = artifactRemoveData(input_data,raw,stimPulseArrs)
    
    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    afterEach(D, @nUpdateWaitbar);
    
    stimPulseArr = stimPulseArrs.singleStimPulseArr;
    
    data_space = size(input_data,1:3);
    nRep = size(input_data,4); nSamples = size(input_data,5);
    chunks = zeros(prod(data_space),nRep,nSamples);
    for idx = 1:prod(data_space)
        [i,j,k] = ind2sub(data_space,idx);
        chunks(idx,:,:) = input_data(i,j,k,:,:); %2d array 
    end
    
    N = prod(data_space);
    p = 1;
    
    detrended_array = zeros(prod(data_space),nRep,nSamples);
    parfor idx = 1:prod(data_space)
        chunk = chunks(idx,:,:);
        detrended_array(idx,:,:) = processChunk(chunk);
        send(D, idx);
    end
    
    copy_raw = raw;
    for idx = 1:prod(data_space)       
        [i,j,k] = ind2sub(data_space,idx);
        temp = detrended_array(idx,:,:);
        for ts = 1:nRep
            stim_times = stimPulseArr{i,k};
            copy_raw(j,stim_times(ts):stim_times(ts)+nSamples-1) = squeeze(temp(1,ts,:));
        end
    end
    
    output_data = copy_raw;
    delete(h);
    
        function nUpdateWaitbar(~)
            waitbar(p/N, h);
            p = p + 1;
        end

end

function dtrend = processChunk(chunk)
    chunk = squeeze(chunk);
    a = chunk - mean(chunk);
    dx = gradient(a,1);
    v = movvar(dx',30);
    cutoff = find(mean(v,2) < 100); 
    assert(isempty(cutoff) == 0,'Threshold for mean windowed variance is too low!')
    cutoff = cutoff(1); 
    blanked = a(:,cutoff:end);
    dtrend = zeros(10,size(a,2));
    for trial = 1:10
        b = blanked(trial,:);
        db = nt_detrend(b',5,[],[],[],[],60)';
        z = 0;
        while ~((db(1) > -35))
            db(1) = [];
            z = z + 1;
        end
        dtrend(trial,:) = [zeros(1,cutoff-1) zeros(1,z) db];
    end
end
% [detrended_array,output_data] = artifactRemoveData(input_data,raw,stimPulseArr,t);
