function varargout = loadData()
%loadData loads extracellular neural recording data and enabled channel indices 

%   loadData() reads an Intan.rhs file and outputs relevant files
if nargout > 2
    % This is the data from the read_Intan_RHS2000 function     
    [raw,amplifier_channels,stim_data,t]=read_Intan_RHS2000_file;
    clear amp_settle_data charge_recovery_data compliance_limit_data...
    spike_triggers frequency_parameters notes reference_channel stim_parameters
    port_number = [amplifier_channels.port_number];
    native_order = [amplifier_channels.native_order]+1;
    enabled_ch_names = {amplifier_channels(:).native_channel_name}';

    % Map individual port indices to 128-order indices using .rhs data
    new_order = [];
    for port_idx = 1:4
        x = port_number == port_idx; 
        new_order  = [new_order native_order(x)+(32*(port_idx-1))];       %#ok<*AGROW>
    end
    varargout = {raw enabled_ch_names new_order stim_data t};
elseif nargout == 2
    % Map using .csv file nad output just enabled channel names and order
    % Useful if you don't want to waste time loading Intan.rhs file
    [dataFile,dataPath] = uigetfile('*.csv','Load .csv file');
    T = readtable(strcat(dataPath,dataFile),'ReadVariableNames',true); 
    enabled = table2array(T(:,4)) == 1;
    enabled_ch_names = table2array(T(enabled,2));

    new_order = [];
    for i = 1:numel(enabled_ch_names)
        port = enabled_ch_names{i}(1);
        native_order = str2double(regexp(enabled_ch_names{i},'\d*','Match')) + 1;
        switch port
            case 'A'
                new_order  = [new_order native_order];       %#ok<*AGROW>
            case 'B'
                new_order  = [new_order native_order+32];       %#ok<*AGROW>
            case 'C'
                new_order  = [new_order native_order+64];       %#ok<*AGROW>
            case 'D'
                new_order  = [new_order native_order+96];       %#ok<*AGROW>
        end
    end
    varargout = {enabled_ch_names new_order};
end

end

