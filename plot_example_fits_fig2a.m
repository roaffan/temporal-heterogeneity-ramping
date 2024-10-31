
%% Set up directories

load_dir = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data\";
load_file = 'alm_dataset5_3150ms.mat';
load(strcat(load_dir,load_file))

% load file with all unit recordings
load(fullfile(load_dir, 'alm_mlfit.mat'),'unit_id','unit_type','time_info',...
    'PSTHs','PSTH_right','PSTH_left','normPSTH',...
    'fSpikes','fModel','mu','tau','sigma','selectivity','Rsquared');

%% Set up colors

color_resolution = 6;
cmp = getPyPlot_cMap('BuPu', color_resolution, [], '"C:\Users\roaffan\Anaconda3\envs\caiman\python.exe"');
color = lighten(cmp(end-1,:),.9);

%% Plot average PSTH and ex-Gaussian fit

time = time_info.tAxisPSTH(find(time_info.tAxisPSTH==0)-3150 : find(time_info.tAxisPSTH==0));
preCue = data.tPreCue;

neurons = [520 397, 227, 146,...
           174, 128, 520 146, 314];

for cell_no = 1:size(neurons,2)
    cell_no = neurons(cell_no);
    
    figure;

    hold on
    plot(time, fSpikes(cell_no,:).*1000, 'Color',lighten([0 0 0],.8), 'LineWidth', .35)
    plot(time, fModel(cell_no,:).*1000, 'Color', color, 'LineStyle', '-', 'LineWidth', 1.75)
    ylabel({'Firing rate'; '(Spikes/s)'})
    format_figure(time);    

    title("Rsqr = "+num2str(Rsquared(cell_no)),'FontName', 'Arial','FontSize',10)

    set(gcf,'color','w')
    set(gcf,'paperunits','inches','units','inches')
    set(gcf,'pos',[3 2.5 1.5486    2.3333])
    set(gcf,'PaperPositionMode','auto')
    fig_dir = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\figures_315s\";
    fig_file = "example_fits_"+num2str(cell_no)+"_psths";
    panel_path = fullfile(fig_dir, fig_file+".pdf");
    print(panel_path, '-dpdf', '-r300')

        
end

%%
function [] = format_figure(T)
    % format the figure
    % add sample, delay and go cue onset
    xlim([min(T)+0.15 max(T)-0.1]);
    yRange = ylim();
    ylim([0, yRange(2)]);
    plot([0         0] ,[0, yRange(2)],'k:','LineWidth',1)
    plot([-3.15 -3.15] ,[0, yRange(2)],'k:','LineWidth',1)
    plot([-2       -2] ,[0, yRange(2)],'k:','LineWidth',1)
    
    xticks([-3.15,-2,0])
    xticklabels({'-3.15','-2','0'})
    yticks([0,yRange(2)/2,round(yRange(2))])
    yticklabels({num2str(0),'',num2str(fix(yRange(2)))})
    
    set(gca, 'FontName', 'Arial','FontSize',9);
    axis square
    box on
end

