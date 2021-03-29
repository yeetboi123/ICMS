function varargout=driver()
    [raw, enabled_ch_names, new_order, stim_data, t, channelInfo, stimInfo, extraInfo]=preprocessData();
    if nargout == 3
        varargout = {channelInfo stimInfo extraInfo};
    else
        varargout = {raw, enabled_ch_names, new_order, stim_data, t, channelInfo, stimInfo, extraInfo};
    end
end

