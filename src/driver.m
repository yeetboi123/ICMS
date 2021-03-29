function varargout=driver()
    [raw, enabled_ch_names, new_order, stim_data, t, channelInfo]=preprocessData();
    if nargout == 1
        varargout{1} = channelInfo;
    end
end