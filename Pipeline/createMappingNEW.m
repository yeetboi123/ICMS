% One mapping function to rule them all...

% Original data presents rows with native channel names ordered by ports
function [raw,stim_data,new_order,stim_new_order,t,save_path] = createMappingNEW
% (A,B,C,D)

% Convert these port order names to 128-channel order names to find that
% channel's location since channel mapping for 128-channel device uses
% 128-channel ordering

% Create map file that maps row to x,y location 

% Produce data structure that contains native channel name, 128-channel
% order, and location data 

%% Assume amplifier_data contains just the channels we are interested in 

[raw,amplifier_channels,stim_data,t]=read_Intan_RHS2000_file;
% This is the data from the read_Intan_RHS2000 function 
port_number = [amplifier_channels.port_number];
native_order = [amplifier_channels.native_order]+1;
native_channel_name = {amplifier_channels(:).native_channel_name}';

% map individual port indices to 128-order indices 
new_order = [];
for port_idx = 1:4
    x = find(port_number == port_idx);   
    new_order  = [new_order native_order(x)+(32*(port_idx-1))];      
end

% map stimulating channels to 128-order indices
stim_chs=find(any(stim_data,2) == 1);
stim_new_order = new_order(stim_chs);
%%  Now find associated location data for each channel using the new order 
save_path = uigetdir('Save path for processing files');

coords = readmatrix('channel_maps/128channel4shank_location_data.csv');

coords_of_interest = coords(new_order,:);

writematrix(coords_of_interest, fullfile(save_path ,'map.csv'));

%% Generate data structure that contains native channel name, 128 channel order, and the channel location 

% I suck with cell arrays so this is the best I could come up with

% First column is the native channel name
% Second column is the 128 channel order
% Third and fourth columns are the x,y coordinates of that channel 
channelInfo = horzcat(native_channel_name, num2cell(new_order'),num2cell(coords_of_interest));
stimChannelInfo = horzcat(native_channel_name(stim_chs), num2cell(stim_new_order'),num2cell(coords_of_interest(stim_chs,:)));

save(fullfile(save_path,'channelInfo'),'channelInfo');
save(fullfile(save_path,'stimChannelInfo'),'stimChannelInfo');
end







