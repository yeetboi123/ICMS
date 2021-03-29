% function plotClusters(processing_path,filtered_data)

addpath mdaio\

input_dir = fullfile(processing_path,'clusterMaps')
mkdir(input_dir);

fs = 3e4;
A = readmda('firings.curated.mda'); % read firings data
primaryCh_vec = A(1,:); 
timeStamp_samples = A(2,:); timeStamp_seconds = timeStamp_samples./fs;
unit_vec = A(3,:);
clear A; 

stimChannelInfo = channelInfo.singleStimCh;
singleRowIdx= stimChannelInfo.singleRowIdx;
single_stim_times = cell2mat(stimInfo(2,:));

nCH = numel(singleRowIdx);

load('channelInfo.mat') % contains info about each channel
allRowIdx=channelInfo.allCh.allRowIdx;
all128OrderIdx=channelInfo.allCh.all128OrderIdx;
allChLocations=channelInfo.allCh.allChLocations;
allChNames=channelInfo.allCh.allChNames;
%%
recording_channels = new_order(allRowIdx);
xlimit1 = -100; xlimit2 = 1000; ylimit1 = -20; ylimit2 = 500;
width = 20; height = 20; spacing = 10; samplesBefore = 15; samplesAfter = 15; scaleFactor = 1;

% xlimit1 = -100; xlimit2 = 1000; ylimit1 = -20; ylimit2 = 500;
% width = 20; height = 20; spacing = 10; samplesBefore = 30; samplesAfter = 30; scaleFactor = 0.2;

for event_idx = 1:200
    ts = timeStamp_samples(event_idx);
    primaryChList = unique(primaryCh_vec); primaryCh = primaryCh_vec(event_idx);
    cluster_idx = unit_vec(event_idx);
    gcf = figure('Menu','none','ToolBar','none');
    set(gcf,'Visible','on')
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
     
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     for ch = primaryChList
%         if ch == primaryCh
%             plot(filtered_data(ch,ts-100:ts+100)','LineWidth', 3);
%         else
%             plot(filtered_data(ch,ts-100:ts+100)','LineWidth', 1); 
%         end
%         hold on;
%     end
%     xline(101)    
    % Plot all channels first

    for i = 1:128
        xlim([xlimit1, xlimit2]);
        ylim([ylimit1, ylimit2]);
        x1 = coords(i,1)-width/2; y1 = coords(i,2)-height/2+spacing;
        x2 = coords(i,1)+width/2; y2 = coords(i,2)+height/2+spacing;
        x = [x1, x2, x2, x1, x1];
        y = [y1, y1, y2, y2, y1];
        xmid = (x1+x2)/2; ymid = (y1+y2)/2;
        condition1 = ismember(i, recording_channels);
        mapped_ch = find(recording_channels == i);
        condition2 = ismember(mapped_ch, primaryChList);

        x_length = xlimit2-xlimit1; y_length = ylimit2-ylimit1;
        if condition1 && condition2
            plot(x, y, '-k','LineWidth', 2);
            hold on
            plot((1:samplesAfter+samplesBefore+1)+x(1)-samplesBefore/2,...
                scaleFactor*filtered_data(mapped_ch,ts-samplesBefore:ts+samplesAfter)'+ymid+height*.25,'LineWidth', 3);
            
        else
            plot(x, y, '-k','LineWidth', 2);
        end

        hold on;
    end
    xlim([xlimit1, xlimit2]);
    ylim([ylimit1, ylimit2]);
    title(sprintf('Cluster %d', cluster_idx))
    uiwait(gcf);
end

%% For a given event, plot the event in a spatial map for each channel

% load('channelInfo.mat') % contains info about each channel
% recording_channels = cell2mat(channelInfo(:,2));
xlimit1 = -100; xlimit2 = 1000; ylimit1 = -20; ylimit2 = 500;
width = 20; height = 20; spacing = 10; samplesBefore = 16; samplesAfter = 16; scaleFactor = 3;
% Plot all channels first

for i = 1:128
    xlim([xlimit1, xlimit2]);
    ylim([ylimit1, ylimit2]);
    x1 = coords(i,1)-width/2; y1 = coords(i,2)-height/2+spacing;
    x2 = coords(i,1)+width/2; y2 = coords(i,2)+height/2+spacing;
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    xmid = (x1+x2)/2; ymid = (y1+y2)/2;
    condition1 = ismember(i, recording_channels);
    mapped_ch = find(recording_channels == i);
    condition2 = ismember(mapped_ch, primaryChList);
    
    x_length = xlimit2-xlimit1; y_length = ylimit2-ylimit1;
    if condition1 && condition2
        plot(x, y, '-k','LineWidth', 2);
        hold on
        plot((1:samplesAfter+samplesBefore+1)+x(1)-samplesBefore/2,...
            scaleFactor*filtered_data(mapped_ch,ts-samplesBefore:ts+samplesAfter)'+ymid + height*.25,'LineWidth', 2);

    else
        plot(x, y, '-k','LineWidth', 2);
    end
    hold on;
end
xlim([xlimit1, xlimit2]);
ylim([ylimit1, ylimit2]);
% If channel in channelInfo, plot the event

% Outline stim channel in red

% Transform coordinates to make it easier to plot waveforms spatially
%%
coords = readmatrix('channel_maps/128channel4shank_location_data.csv');
coords_of_interest = coords(new_order,:);
writematrix(coords_of_interest, fullfile(processing_path ,'map.csv'));

coords_scale = 1
xlimit1 = -100; xlimit2 = 1000; ylimit1 = -20; ylimit2 = 500;
width = 70; height = 20; spacing = 10; samplesBefore = 15; samplesAfter = 15; scaleFactor = 2;
intrashank_space = 80;
for i = 1:128
    xlim([xlimit1, xlimit2]);
    ylim([ylimit1, ylimit2]);
    x1 = coords_scale*coords(i,1)-width/2; y1 = coords_scale*coords(i,2)-height/2+spacing;
    x2 = coords_scale*coords(i,1)+width/2; y2 = coords_scale*coords(i,2)+height/2+spacing;
    [x,xmid] = coord_transform(coords(i,1),width,intrashank_space)
    y = [y1, y1, y2, y2, y1];
%     xmid = (x1+x2)/2; 
    ymid = (y1+y2)/2;
    condition1 = ismember(i, recording_channels);
    mapped_ch = find(recording_channels == i);
    condition2 = ismember(mapped_ch, primaryChList);
    
    x_length = xlimit2-xlimit1; y_length = ylimit2-ylimit1;
    if condition1 && condition2
        plot(x, y, '-k','LineWidth', 2);
        hold on
        plot((1:samplesAfter+samplesBefore+1)+x(1),...
            scaleFactor*filtered_data(mapped_ch,ts-samplesBefore:ts+samplesAfter)'+ymid,'LineWidth', 3);

    else
        plot(x, y, '-k','LineWidth', 2);
    end
    hold on;
end
xlim([xlimit1, xlimit2]);
ylim([ylimit1, ylimit2]);
% end
%%
function [x,xmid] = coord_transform(x_old,width,intrashank_space)
    if (x_old == 30) || (x_old == 330) || (x_old == 630) || (x_old == 930)
        temp = (x_old - 30) + intrashank_space;
    else
        temp = x_old;
    end
    x1 = temp-width/2; 
    x2 = temp+width/2; 
    xmid = (x1+x2)/2
    x = [x1, x2, x2, x1, x1];
end