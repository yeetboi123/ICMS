% Stimulating channels digraph 
function generateStimulationDigraph(stim_graph,nSCH,nCU,singleRowIdx,coordinates,stimChannelInfo,processing_path)

stimDigraphPath = fullfile(processing_path,'stimulationDigraphs');
if ~exist(stimDigraphPath, 'dir')
    mkdir(stimDigraphPath)
    addpath(stimDigraphPath)
end
      

% iterate through each row of the stim graph
% each row corresponds to a stimulating channel's activations on all
% channels
good_stim_ch_num = 0; good_idx_arr = [];
for ii = 1:nSCH
    if ~all(stim_graph(ii,:)==0)
        good_stim_ch_num = good_stim_ch_num + 1;
        good_idx_arr = [good_idx_arr ii];
    end
end

colors=vals2colormap([1:nCU],'lines');
pos = cell2mat(coordinates);
stim_ch_names = stimChannelInfo.singleChNames;

counter = good_stim_ch_num;
new_stim_graph = zeros(size(stim_graph));
while counter > 0
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:min(9,counter)

        subplot(3,3,i)
        good_idx = good_idx_arr(i);
        new_stim_graph(good_idx,:) = stim_graph(good_idx,:); % insert row of interest in zero matrix
        G = digraph(new_stim_graph); G = flipedge(G); % flips edge
        h=plot(G,'Xdata',pos(:,1),'YData',pos(:,2),'ArrowSize',10);
        
        stimulating_ch = singleRowIdx(good_idx);
        highlight(h,stimulating_ch,'NodeColor','r')

        % color edges based on edge weights
        custom_green = [64, 226, 15]/255;
        for edge = 1:numel(G.Edges.Weight)
            highlight(h,G.Edges.EndNodes(edge,:),'EdgeColor',colors(G.Edges.Weight(edge),:))
            highlight(h,G.Edges.EndNodes(edge,1),'NodeColor',custom_green)
            hold on;
        end
        
        name = stim_ch_names{good_idx};
        h.LineWidth=2;
        h.Marker='s';
        h.MarkerSize=4;
        h.NodeLabelMode = 'auto'; 
        h.NodeLabel = {};
        title(sprintf('Minimal current for activation digraph for stim channel %s',name))

        % generate dummy date for legend 
        hold on;
        l1 = scatter(nan, nan, [],colors(1,:),'filled','^');
        l2 = scatter(nan, nan, [],colors(2,:),'filled','^');
        l3 = scatter(nan, nan, [],colors(3,:),'filled','^');
        l4 = scatter(nan, nan, [],colors(4,:),'filled','^');
        lgd = legend([l1, l2, l3, l4], {'0.5 uA', '1 uA', '2 uA', '5 uA'}); % hardcoded!
        lgd.FontSize = 7;
        lgd.Location = 'southwest';

        ylim([-140 500]); % ensures legend is not covering contacts

        % get min and max x and y coordinates of each shank (x1,x2,y1,y2)
        shifts = [-30 -18 60 36]; % shift_array
        shank1_border = [0 0 30 450]; rectangle('Position',shank1_border + shifts,'Linewidth',2,'Edgecolor',colors(1,:))
        shank2_border = [300 0 30 450]; rectangle('Position',shank2_border + shifts,'Linewidth',2,'Edgecolor',colors(1,:))
        shank3_border = [600 0 30 450]; rectangle('Position',shank3_border + shifts,'Linewidth',2,'Edgecolor',colors(1,:))
        shank4_border = [900 0 30 450]; rectangle('Position',shank4_border + shifts,'Linewidth',2,'Edgecolor',colors(1,:))

        set(gca,'xtick',[]); set(gca,'ytick',[])
    end
    
    exportgraphics(gcf, 'Stimulation digraphs for ICMS-15 Feb.13')

    counter = counter - 9;
end

end