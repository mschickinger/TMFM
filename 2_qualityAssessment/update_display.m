
%% red
ch = 1;
%r/rms
spRMS = subplot('Position', [trace_left 2*h0+0.1*h0+5*trace_height trace_width trace_height]);
tmp_plot = (1:size(data{m}{s,ch}.vwcm.pos,1))';
tmp_plot = [tmp_plot(data{m}{s,ch}.vwcm.pos(:,1)>0) data{m}{s,ch}.vwcm.r(data{m}{s,ch}.vwcm.pos(:,1)>0)];
hold off
if any(tmp_plot(:,2)>5)
    for j = find(tmp_plot(:,2)'>5)'
        plot(tmp_plot(j,1)*[1 1], [0 pxLIM], '-', 'Color', [1 1 1].*.7);
        hold on
    end
end
plot(tmp_plot(:,1), tmp_plot(:,2),['.' Colors{ch}], 'MarkerSize', 8);
hold on
plot(tmp_plot(:,1),data{m}{s,ch}.vwcm.rms10(data{m}{s,ch}.vwcm.pos(:,1)>0),'-k', 'LineWidth', 0.5);

ax = gca;
ax.XLim = [0 XT(end)];
ax.YLim = [0 pxLIM];
ax.XTick = XT;
ax.XTickLabel = {''};
ax.Layer = 'top';
grid on
box off
title(['Movie ' num2str(m) ', spot pair ' num2str(s) ' of ' num2str(size(data{m},1)) '.'], 'FontSize', 12, 'Interpreter', 'tex')

% red intensity
spINT = subplot('Position', [trace_left 2*h0+4*trace_height trace_width trace_height]);
hold off
plot(data{m}{s,ch}.itrace,'Color',Colors{ch},'LineWidth',1)
hold on
plot(data{m}{s,ch}.med_itrace,'k','LineWidth',.5)
plot(data{m}{s,ch}.med_itrace,'k','LineWidth',.5)
plot(fit_cutoff{m,ch}(s)*[1 1],[0 max(data{m}{s,ch}.itrace)],'k--')
ax = gca;
ax.XLim = [0 XT(end)];
ax.XTick = XT;
ax.XTickLabel = XTL;
ax.YLim = [0.95*min(data{m}{s,ch}.itrace) max(data{m}{s,ch}.itrace)];
if str2double(ax.YTickLabel{1}) == 0.95*min(data{m}{s,ch}.itrace)
    ax.YTickLabel{1} = '';
end
if str2double(ax.YTickLabel{end}) == max(data{m}{s,ch}.itrace)
    ax.YTickLabel{end} = '';
end
ax.Layer = 'top';
grid on
box off

% red scatter
spSC = subplot('Position', [1.5*trace_left+trace_width spRMS.Position(2) trace_height*scrsz(4)/scrsz(3) trace_height]);
hold off
plot(0, 0, [Colors{ch} 'x'])
plot(data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r>5,1),data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r>5,2), 'o', 'Color', [1 1 1]*.7, 'MarkerSize', 3);
hold on
plot(data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r<=5,1),data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r<=5,2), [Colors{ch} '+'], 'MarkerSize', 3);
scLIM = max(5,2*ceil((pxLIM+2)/2));
xlim([-scLIM scLIM])
ylim([-scLIM scLIM])
TC = cell(1,2*scLIM+1); TC{1} = -scLIM; TC{(length(TC)+1)/2} = 0; TC{end} = scLIM;
set(gca, 'XTick', -scLIM:scLIM, 'YTick', -scLIM:scLIM, 'XTickLabel', TC, 'YTickLabel', TC);
grid on, axis square
title('Positions (drift corrected)')

% red pos/pos0
subplot('Position', [sum(spSC.Position([1 3 3]))+trace_left spINT.Position(2) 1-(sum(spSC.Position([1 3 3]))+1.5*trace_left) 1-spINT.Position(2)-h0]);
hold off
plot(data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),1),data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),2),[Colors{ch} '.'],'MarkerSize',8)
hold on
plot(data{m}{s,ch}.vwcm.medians101(:,1),data{m}{s,ch}.vwcm.medians101(:,2),'Color', [1 1 1]*.7,'LineWidth', 2)
plot(data{m}{s,ch}.pos0(tmp_plot(:,1),1),data{m}{s,ch}.pos0(tmp_plot(:,1),2),'k-x','MarkerSize', 16, 'LineWidth', 2)
axis equal
if size(data{m}{s,ch}.pos0(tmp_plot(:,1),:),1)>0
    xlim([min(data{m}{s,ch}.pos0(tmp_plot(:,1),1))-5 max(data{m}{s,ch}.pos0(tmp_plot(:,1),1))+5])
    ylim([min(data{m}{s,ch}.pos0(tmp_plot(:,1),2))-5 max(data{m}{s,ch}.pos0(tmp_plot(:,1),2))+5])
end
grid on
set(gca, 'YDir', 'normal')
title('Pos and Pos0')

% red deviations pos from pos0
posSC = subplot('Position', [sum(spSC.Position([1 3]))+0.5*trace_left spRMS.Position(2) spSC.Position(3:4) ]);
hold off
plot(data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),1)-data{m}{s,ch}.pos0(tmp_plot(:,1),1),data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),2)-data{m}{s,ch}.pos0(tmp_plot(:,1),2),[Colors{ch} '.'],'MarkerSize',8)
axis equal
xlim([-scLIM scLIM])
ylim([-scLIM scLIM])
grid on
set(gca, 'YDir', 'normal', 'XTick', -scLIM:scLIM, 'YTick', -scLIM:scLIM, 'XTickLabel', TC, 'YTickLabel', TC);          
title('Pos - Pos0')

% red start/end pos avg
% start
subplot('Position',[spSC.Position(1), spINT.Position(2)-.02, trace_height*scrsz(4)/scrsz(3), trace_height])
hold off
%x_0 = round(mean(data{m}{s,ch}.pos0(1:100,1)));
%y_0 = round(mean(data{m}{s,ch}.pos0(1:100,2)));
if size(avg_img,2) == 4
    plot_subframe(avg_img{m,ch}, x0{m}(s,1), y0{m}(s,1), 12);
elseif size(avg_img,2) == 2
    plot_subframe(avg_img{m,ch}{1}, x0{m}(s,1), y0{m}(s,1), 12);
end
hold on        
ellipse(r_integrate, r_integrate, 0, x0{m}(s,1), y0{m}(s,1), Colors{ch});
title(['Start Pos: ' num2str(x0{m}(s,1)) ', ' num2str(y0{m}(s,1))]);
axis square
set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});

% end
subplot('Position',[posSC.Position(1), spINT.Position(2)-.02, trace_height*scrsz(4)/scrsz(3), trace_height])
hold off
%x_0 = round(mean(data{m}{s,ch}.pos0(end-99:end,1)));
%y_0 = round(mean(data{m}{s,ch}.pos0(end-99:end,2)));
if size(avg_img,2) == 4
    plot_subframe(avg_img{m,ch+2}, x0{m}(s,2), y0{m}(s,2), 12);
elseif size(avg_img,2) == 2
    plot_subframe(avg_img{m,ch}{end}, x0{m}(s,2), y0{m}(s,2), 12);
end
hold on        
ellipse(r_integrate, r_integrate, 0, x0{m}(s,2), y0{m}(s,2), Colors{ch});
title(['End Pos: ' num2str(x0{m}(s,2)) ', ' num2str(y0{m}(s,2))]);
axis square
set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});

%% green 
ch = 2;
% green r/rms
spRMS = subplot('Position', [trace_left h0+0.3*h0+3*trace_height trace_width trace_height]);
tmp_plot = (1:size(data{m}{s,ch}.vwcm.pos,1))';
tmp_plot = [tmp_plot(data{m}{s,ch}.vwcm.pos(:,1)>0) data{m}{s,ch}.vwcm.r(data{m}{s,ch}.vwcm.pos(:,1)>0)];
hold off
if any(tmp_plot(:,2)>5)
    for j = find(tmp_plot(:,2)>5)'
        plot(tmp_plot(j,1)*[1 1], [0 pxLIM], '-', 'Color', [1 1 1].*.7);
        hold on
    end
end
plot(tmp_plot(:,1), tmp_plot(:,2),['.' Colors{ch}], 'MarkerSize', 8);
hold on
plot(tmp_plot(:,1),data{m}{s,ch}.vwcm.rms10(data{m}{s,ch}.vwcm.pos(:,1)>0),'-k', 'LineWidth', 0.5);

ax = gca;
ax.XLim = [0 XT(end)];
ax.YLim = [0 pxLIM];
ax.XTick = XT;
ax.XTickLabel = {''};
ax.Layer = 'top';
grid on
box off

% green intensity
spINT = subplot('Position', [trace_left h0+0.2*h0+2*trace_height trace_width trace_height]);
hold off
plot(data{m}{s,ch}.itrace,'Color',Colors{ch},'LineWidth',1)
hold on
plot(data{m}{s,ch}.med_itrace,'k','LineWidth',.5)
plot(fit_cutoff{m,ch}(s)*[1 1],[0 max(data{m}{s,ch}.itrace)],'k--')
ax = gca;
ax.XLim = [0 XT(end)];
ax.XTick = XT;
ax.XTickLabel = {''};
ax.YLim = [0.95*min(data{m}{s,ch}.itrace) max(data{m}{s,ch}.itrace)];
if str2double(ax.YTickLabel{1}) == 0.95*min(data{m}{s,ch}.itrace)
    ax.YTickLabel{1} = '';
end
if str2double(ax.YTickLabel{end}) == max(data{m}{s,ch}.itrace)
    ax.YTickLabel{end} = '';
end
ax.Layer = 'top';
grid on
box off

% green scatter
spSC = subplot('Position', [1.5*trace_left+trace_width spRMS.Position(2) trace_height*scrsz(4)/scrsz(3) trace_height]);
hold off
plot(0, 0, [Colors{ch} 'x'])
plot(data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r>5,1),data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r>5,2), 'o', 'Color', [1 1 1]*.7, 'MarkerSize', 3);
hold on
plot(data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r<=5,1),data{m}{s,ch}.vwcm.dispmed101(data{m}{s,ch}.vwcm.r<=5,2), [Colors{ch} '+'], 'MarkerSize', 3);
%scLIM = max(5,pxLIM+2);
xlim([-scLIM scLIM])
ylim([-scLIM scLIM])
TC = cell(1,2*scLIM+1); TC{1} = -scLIM; TC{(length(TC)+1)/2} = 0; TC{end} = scLIM;
set(gca, 'XTick', -scLIM:scLIM, 'YTick', -scLIM:scLIM, 'XTickLabel', TC, 'YTickLabel', TC);
grid on, axis square
title('Positions (drift corrected)')

% green pos/pos0

subplot('Position', [sum(spSC.Position([1 3 3]))+trace_left spINT.Position(2) 1-(sum(spSC.Position([1 3 3]))+1.5*trace_left) sum(spRMS.Position([1 3]))-spINT.Position(2)-h0]);
hold off
plot(data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),1),data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),2),[Colors{ch} '.'],'MarkerSize',8)
hold on
plot(data{m}{s,ch}.vwcm.medians101(:,1),data{m}{s,ch}.vwcm.medians101(:,2),'Color', [1 1 1]*.7,'LineWidth', 2)
plot(data{m}{s,ch}.pos0(tmp_plot(:,1),1),data{m}{s,ch}.pos0(tmp_plot(:,1),2),'k-x','MarkerSize', 16, 'LineWidth', 2)
axis equal
if size(data{m}{s,ch}.pos0(tmp_plot(:,1),:),1)>0
    xlim([min(data{m}{s,ch}.pos0(tmp_plot(:,1),1))-5 max(data{m}{s,ch}.pos0(tmp_plot(:,1),1))+5])
    ylim([min(data{m}{s,ch}.pos0(tmp_plot(:,1),2))-5 max(data{m}{s,ch}.pos0(tmp_plot(:,1),2))+5])
end
grid on
set(gca, 'YDir', 'normal')

% green deviations pos from pos0
posSC = subplot('Position', [sum(spSC.Position([1 3]))+0.5*trace_left spRMS.Position(2) spSC.Position(3:4) ]);
hold off
plot(data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),1)-data{m}{s,ch}.pos0(tmp_plot(:,1),1),data{m}{s,ch}.vwcm.pos(tmp_plot(:,1),2)-data{m}{s,ch}.pos0(tmp_plot(:,1),2),[Colors{ch} '.'],'MarkerSize',8)
axis equal
xlim([-scLIM scLIM])
ylim([-scLIM scLIM])
grid on
set(gca, 'YDir', 'normal', 'XTick', -scLIM:scLIM, 'YTick', -scLIM:scLIM, 'XTickLabel', TC, 'YTickLabel', TC);          
title('Pos - Pos0')

% green start/end pos avg
% start
subplot('Position',[spSC.Position(1), spINT.Position(2)-.02, trace_height*scrsz(4)/scrsz(3), trace_height])
hold off
%x_0 = round(mean(data{m}{s,ch}.pos0(1:100,1)));
%y_0 = round(mean(data{m}{s,ch}.pos0(1:100,2)));
if size(avg_img,2) == 4
    plot_subframe(avg_img{m,ch}, x0{m}(s,3), y0{m}(s,3), 12);
elseif size(avg_img,2) == 2
    plot_subframe(avg_img{m,ch}{1}, x0{m}(s,3), y0{m}(s,3), 12);
end
hold on        
ellipse(r_integrate, r_integrate, 0, x0{m}(s,3), y0{m}(s,3), Colors{ch});
title(['Start Pos: ' num2str(x0{m}(s,3)) ', ' num2str(y0{m}(s,3))]);
axis square
set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});

% end
subplot('Position',[posSC.Position(1), spINT.Position(2)-.02, trace_height*scrsz(4)/scrsz(3), trace_height])
hold off
%x_0 = round(mean(data{m}{s,ch}.pos0(end-99:end,1)));
%y_0 = round(mean(data{m}{s,ch}.pos0(end-99:end,2)));
if size(avg_img,2) == 4
    plot_subframe(avg_img{m,ch+2}, x0{m}(s,4), y0{m}(s,4), 12);
elseif size(avg_img,2) == 2
    plot_subframe(avg_img{m,ch}{end}, x0{m}(s,4), y0{m}(s,4), 12);
end
hold on        
ellipse(r_integrate, r_integrate, 0, x0{m}(s,4), y0{m}(s,4), Colors{ch});
title(['End Pos: ' num2str(x0{m}(s,4)) ', ' num2str(y0{m}(s,4))]);
axis square
set(gca, 'YDir', 'normal', 'XTickLabel', {}, 'YTickLabel', {});

%% displacements
%y
subplot('Position', [trace_left h0 trace_width trace_height])
hold off
plot(disp_from_map{m}{s}(:,2), '.', 'MarkerSize', 8)
hold on
plot(disp_median_from_map{m}{s}(:,2), 'Color', 'c', 'LineWidth', .5)
ax = gca;
ax.XTick = XT;
ax.XTickLabel = XTL;
ax.YLim = [-pxLIM pxLIM];
ax.YTick = -pxLIM:pxLIM;
ax.Layer = 'top';
ylabel('deltaY/medY')
grid on
box off

%x
subplot('Position', [trace_left h0+0.1*h0+trace_height trace_width trace_height])
hold off
plot(disp_from_map{m}{s}(:,1), '.', 'MarkerSize', 8)
hold on
plot(disp_median_from_map{m}{s}(:,1), 'Color', 'c', 'LineWidth', .5)
ax = gca;
ax.XTick = XT;
ax.XTickLabel = {};
ax.YLim = [-pxLIM pxLIM];
ax.YTick = -pxLIM:pxLIM;
ax.Layer = 'top';
ylabel('deltaX/medX')
grid on
box off

%scatter
subplot('Position', [spSC.Position(1)+0.01 h0 bg.Position(1)-spSC.Position(1)-0.02 spINT.Position(2)-h0-0.02])
plot(disp_from_map{m}{s}(:,1),disp_from_map{m}{s}(:,2),'.','MarkerSize',8)
axis square
ax = gca;
ax.XLim = [-scLIM scLIM];
ax.YLim = [-scLIM scLIM];
set(gca, 'YDir', 'normal', 'XTick', -scLIM:scLIM, 'YTick', -scLIM:scLIM, 'XTickLabel', TC, 'YTickLabel', TC);
ax.Layer = 'top';
grid on
box off