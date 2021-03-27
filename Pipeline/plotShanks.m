% Grab enabled channels only from impedance file
addpath channel_maps/
[dataFile,dataPath] = uigetfile('*.csv');
T = readtable(strcat(dataPath,dataFile),'ReadVariableNames',true); 
channel_names = table2array(T(:,1));
enabled = find(table2array(T(:,4)));

new_order = enabled
impedances = table2array(T(:,5));

coords=readmatrix(strcat('channel_maps/128channel4shank_location_data.csv'));
shank1 = [17:32 81:96];
shank2 = [1:16 65:80];
shank3 = [33:48 97:112];
shank4 = [49:64 113:128];
shank_combined = {shank1 shank2 shank3 shank4};
xlimit1 = -100; xlimit2 = 1000; ylimit1 = -20; ylimit2 = 500;
width = 20; height = 20; spacing = 10

[Lia1,Locb1]=ismember(shank1,new_order);
[Lia2,Locb2]=ismember(shank2,new_order);
[Lia3,Locb3]=ismember(shank3,new_order);
[Lia4,Locb4]=ismember(shank4,new_order);
arr = {Locb1,Locb2,Locb3,Locb4};
for shank = 1:4
    shank_ind = shank_combined{shank};
    shank_coords = coords(shank_ind,:);
    stim_ind = []; stim_ch_flag = 0;

    shank_arr = arr{shank}
    
    for i = 1:32
        
        if shank_arr(i) > 0 && ismember(new_order(shank_arr(i)),shank_ind)
            stim_ch_flag = 1; 
        else
            stim_ch_flag = 0;
        end
        
        x1 = shank_coords(i,1)-width/2; y1 = shank_coords(i,2)-height/2+spacing;
        x2 = shank_coords(i,1)+width/2; y2 = shank_coords(i,2)+height/2+spacing;

        x = [x1, x2, x2, x1, x1];
        y = [y1, y1, y2, y2, y1];
        xmid = (x1+x2)/2; ymid = (y1+y2)/2;
        if stim_ch_flag
            plot(x, y, '-r','LineWidth', 2);
            text(x1+1, ymid, channel_names(new_order(shank_arr(i))),'FontSize',7)
%             text(x1, ymid+15, num2str(impedances(new_order(shank_arr(i)))),'FontSize',7)
            
        else
            plot(x, y, '-k','LineWidth', 2);
        end
        hold on;

        xlim([xlimit1, xlimit2]);
        ylim([ylimit1, ylimit2]);
        
    end
    hold on;
end
title('Mar 9: ICMS-15 spatial layout of stimulation channels')
