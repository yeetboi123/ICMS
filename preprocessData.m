function [raw enabled_ch_names new_order stim_data t]=preprocessData()
    
    addpath util\
    [raw enabled_ch_names new_order stim_data t]=loadData();
    clear amp_settle_data charge_recovery_data compliance_limit_data...
    spike_triggers frequency_parameters notes reference_channel stim_parameters

%     [enabled_ch_names new_order]=loadData;
    processing_path = uigetdir('Select directory to save intermediary files for preprocessing');
    


end