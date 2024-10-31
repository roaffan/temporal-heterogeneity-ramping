
clear
plot = 0;

%% Load data

dir_name = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data\FixedDelayTask\";
files = dir(fullfile(dir_name, '*.mat'));
myTable = readtable(dir_name + "SI_table_2_bialteral_perturb.xlsx",...
              'PreserveVariableNames',true);
          
units_N = double(table2array(myTable([1:size(files,1)]+1,'All units')));

%% Initialize varaiables 

psth_right = []; 
psth_left  = [];
time_axis  = [];

sample_onset = [];
delay_onset  = [];
go_cue_onset = [];

spike_rate_pre   = [];
spike_rate_delay = [];
spike_rate_post  = [];
spikes = [];

cell_type = [];
trial_type = [];
photo_stim = [];

delay_activity = [];
delay_selectivity = [];
norm_delay_selectivity = [];

%% Get PSTHs for correct right vs left lick trials

unit_j = 0;
% Loop through files
for file_i = 1:size(files,1)
    file_name = fullfile(dir_name, files(file_i).name);
    load(file_name)
    N = units_N(file_i);
    % Loop through units within each file
    for unit_i = 1:N
        j = unit_j + unit_i;

        [right, left, T, unitInfo] = get_psth_all_755_units(unit(unit_i));
        psth_right = [psth_right; right];
        psth_left  = [psth_left; left];
   
        time_axis    = [time_axis; T.tAxisPSTH];
        sample_onset = [sample_onset; T.meanSampleOnset - T.meanCueOnset];
        delay_onset  = [delay_onset; T.meanDelayOnset - T.meanCueOnset];
        go_cue_onset = [go_cue_onset; T.meanCueOnset];
        
        profile(j).cell_type        = [unitInfo.type];
        profile(j).depth            = [unitInfo.depth];

        profile(j).spike_rate_pre   = [spike_rate_pre; unitInfo.spikes_presample];
        profile(j).spike_rate_delay = [spike_rate_delay; unitInfo.spikes_delay];
        profile(j).spike_rate_post  = [spike_rate_post; unitInfo.spikes_postcue];

        profile(j).spike_rate_pre   = [spike_rate_pre; unitInfo.spikes_presample];
        profile(j).spike_rate_delay = [spike_rate_delay; unitInfo.spikes_delay];
        profile(j).spike_rate_post  = [spike_rate_post; unitInfo.spikes_postcue];
        profile(j).spikes           = [spikes; unitInfo.spikes'];

        profile(j).trial_type       = [trial_type; unitInfo.trial_type'];
        profile(j).photo_stim       = [photo_stim; unitInfo.photo_stim'];
        
        profile(j).delay_selective        = unitInfo.selective;
        profile(j).delay_preference       = unitInfo.delay_preference;
        profile(j).norm_delay_selectivity = unitInfo.norm_delay_selectivity;

        profile(j).delay_sig      = unitInfo.delay_sig;
        profile(j).ramp_direction = unitInfo.ramp_direction;

        profile(j).PSTH     = unitInfo.PSTH;
        profile(j).normPSTH = unitInfo.normPSTH;
    end
    unit_j = unit_j+N;
end

%% Save data
time_info = T;

filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";
save(fullfile(filename, 'unit_profile.mat'),'profile','psth_left',...
    'psth_right','time_info','-v7.3');

%% Get index of unit types

RS_ind = [profile.cell_type]=="RS";
FS_ind = [profile.cell_type]=="FS";
DA_ind = [profile.delay_sig];
up_ind = [profile.ramp_direction]==1;
down_ind = [profile.ramp_direction]==-1;

selective_ind = [profile.delay_selective];
rigt_ind = [profile.norm_delay_selectivity] > 0;
left_ind = [profile.norm_delay_selectivity] < 0;

delay_active_RS = sum(RS_ind & DA_ind);
delay_active_FS = sum(FS_ind & DA_ind);

up_RS = sum((RS_ind & DA_ind) & up_ind);
up_FS = sum(FS_ind & DA_ind & up_ind);

down_RS = sum(RS_ind & DA_ind & down_ind);
down_FS = sum(FS_ind & DA_ind & down_ind);

index.pyramidal  = [RS_ind];
index.fastspike  = [FS_ind];
index.ramp_up    = [DA_ind & up_ind];
index.ramp_down  = [DA_ind & down_ind];
index.trial_type = [selective_ind];

filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";
save(fullfile(filename, 'alm_index.mat'),'index','-v7.3');

%% Report counts of ramp up and ramp down putative RS and FS neurons

disp(['Number of RS active during delay: ', num2str(delay_active_RS), ' units, ',...
      num2str((delay_active_RS/667)*100), '% of 667 units']);
disp(['Number of FS active during delay: ', num2str(delay_active_FS), ' units, ',...
      num2str((delay_active_FS/74)*100), '% of 74 units']);

disp(['Number of RS ramp up: ', num2str(up_RS), ' units, ',...
      num2str((up_RS/delay_active_RS)*100), '% of ', num2str(delay_active_RS),...
      ' units']);
disp(['Number of FS ramp up: ', num2str(up_FS), ' units, ',...
      num2str((up_FS/delay_active_FS)*100), '% of ', num2str(delay_active_FS),...
      ' units']);

disp(['Number of RS ramp down: ', num2str(down_RS), ' units, ',...
      num2str((down_RS/delay_active_RS)*100), '% of ', num2str(delay_active_RS),...
      ' units']);
disp(['Number of FS ramp down: ', num2str(down_FS), ' units, ',...
      num2str((down_FS/delay_active_FS)*100), '% of ', num2str(delay_active_FS),...
      ' units']);

%% Plot PSTHs of specific unit types

if plot

go_cue_time = find(time_info.tAxisPSTH(1,:)==0);
trial_onset_time = find(time_info.tAxisPSTH(1,:)==0) - 3150;

figure;

    for k = 1:4
        switch k
            case 1
                UOI_index = find(RS_ind & DA_ind & up_ind);
            case 2
                UOI_index = find(FS_ind & DA_ind & up_ind);
            case 3
                UOI_index = find(RS_ind & DA_ind & down_ind);
            case 4
                UOI_index = find(FS_ind & DA_ind & down_ind);
        end
    
    x = time_info.tAxisPSTH(:, trial_onset_time : go_cue_time);
    x = repmat(x,[755,1]);
    psth = (psth_right(:, trial_onset_time : go_cue_time)+...
            psth_left(:, trial_onset_time : go_cue_time))./ 2;
    y = normpsth;
    
    for i = 1:length(UOI_index)
        subplot(2,2,k)
        hold on
        plot(x(UOI_index(i),:), y(UOI_index(i),:))
    %     plot(x(20,:), y(20,:))
    end
    
    % xlim([time_info.meanSampleOnset- time_info.meanCueOnset time_info.meanCueOnset- time_info.meanCueOnset]);
    yRange = ylim();
    plot([0 0],[-1000 1000],'k:',"LineWidth",1)
    plot([time_info.meanSampleOnset - time_info.meanCueOnset  time_info.meanSampleOnset - time_info.meanCueOnset] ,yRange,'k:',"LineWidth",1)
    plot([time_info.meanDelayOnset - time_info.meanCueOnset    time_info.meanDelayOnset - time_info.meanCueOnset]   ,yRange,'k:',"LineWidth",1)
    plot([-1.5750    -1.5750]   ,yRange,'k--',"LineWidth",1)
    xlabel('Time from go cue osnet (s)')
    % ylabel('(SR delay - avg SR baseline) / avg SR baseline')
    ylabel('\Delta norm. Firing Rate [Hz]')
    ylim(yRange);
    
    switch k
        case 1
            title("Ramp-up pyramidal n = 277")
        case 2
            title("Ramp-up fast-spiking n = 48")
        case 3
            title("Ramp-down pyramidal n = 173")
        case 4
            title("Ramp-down fast-spiking n = 19")
    end
    
    % axis square
    set(gcf, 'Position',  [300   116   859   717])
    end

end

%% Function modified from original code by Inagaki et al., 2018 (retrieved from CRCNS)

function [psth_correct_right, psth_correct_left, T, unit_info ] = get_psth_all_755_units(unit)

    % set parameters for plot
    preCueDur  = 4.2; % in sec: plot how many seconds before go cue onset. 
    postCueDur = 1.5; % in sec: plot how many seconds after go cue onset. 
    
    trialType = unit.Behavior.Trial_types_of_response_vector; % trial type based on outcome
    % 1: correct lick R trial,  2: correct lick L trial,  3: incorrect lick R trial,  4: incorrect lick L trial,
    % 5: correct lick R trial with early lick,   6: correct lick L trial with early lick, 
    % 7: incorrect lick R trial with early lick, 8: incorrect lick L trial with early lick
    % 9: lick R trial but no response,          10: lick L trial but no response,
    % 11: others(unidentified)
    
    photoStimTrial = unit.Behavior.stim_trial_vector;  % photo stim trial type
    % 0: no stim,   1: 0.05mW bilateral stim during early dealy,  
    % 2: 0.1mW bilateral stim during early dealy,
    % 3: 0.2mW bilateral stim during early dealy,
    % 4: 0.3mW bilateral stim during early dealy,
    % don't analyze 0.05mW it was too weak to have any behavioral effect
    
    % extract timing info
    T.tAxisPSTH   = -preCueDur-0.1:0.001:postCueDur+0.1; % T axis for PSTH, shift by 0.1 to remove smoothing artifact at the edge
    T.sampleOnset = [unit.Behavior.Sample_start]; % smaple epoch onset
    T.delayOnset  = [unit.Behavior.Delay_start]; % delay epoch onset
    T.cueOnset    = [unit.Behavior.Cue_start]; % go cue  onset
     
    % calculate mean onset time of each epoch for plotting   
    % exclude early lick and non-response trial by "trialType<5"
    T.meanSampleOnset = mean(T.sampleOnset(trialType<5));
    T.meanDelayOnset  = mean(T.delayOnset(trialType<5));
    T.meanCueOnset    = mean(T.cueOnset(trialType<5));
    
    % extract spikes  from each trial
    spikeTimes       = unit.SpikeTimes; % time bin of spike
    trialIdxOfSpikes = unit.Trial_idx_of_spike; % trial num of spike
    
    trialRange            = unit.Trial_info.Trial_range_to_analyze; % range of trials to analyze
    trialTypeInRange      = trialType(trialRange(1) : trialRange(2));        % trial type in range
    photoStimTrialInRange = photoStimTrial(trialRange(1) : trialRange(2));   % stim trial type in range
    
    PSTH         = nan( trialRange(2)- trialRange(1)+1,numel(T.tAxisPSTH));
    SR_presample = nan( trialRange(2)- trialRange(1)+1, 1);
    SR_sample    = nan( trialRange(2)- trialRange(1)+1, 1);
    SR_delay     = nan( trialRange(2)- trialRange(1)+1, 1);
    SR_postcue   = nan( trialRange(2)- trialRange(1)+1, 1);
    
    for tr = trialRange(1) : trialRange(2)
    
        % extarct spikes of each trial
        SpikesTmp   = spikeTimes(trialIdxOfSpikes == tr) - T.cueOnset(tr); % align to go cue 
        PSTHTmp     = hist(SpikesTmp,T.tAxisPSTH)*1000;
        PSTH(tr - trialRange(1) +1,:)  = smooth(PSTHTmp,101);
                
        sample_onset = T.sampleOnset(tr - trialRange(1) +1) - T.cueOnset(tr);
        delay_onset  = T.delayOnset(tr - trialRange(1) +1) - T.cueOnset(tr);
        response_onset = T.cueOnset(tr - trialRange(1) +1) - T.cueOnset(tr);
    
        SR_presample(tr - trialRange(1) +1,:) = sum(SpikesTmp < sample_onset) ./ (preCueDur - abs(sample_onset));
        SR_sample(tr - trialRange(1) +1,:)    = sum(SpikesTmp > sample_onset & SpikesTmp < delay_onset) ./ (delay_onset-sample_onset);
        SR_delay(tr - trialRange(1) +1,:)     = sum(SpikesTmp > delay_onset & SpikesTmp < response_onset) ./ (response_onset-delay_onset);
        SR_delay1(tr - trialRange(1) +1,:)    = sum(SpikesTmp > delay_onset & SpikesTmp < response_onset-1) ./ ((response_onset-1)-delay_onset);
        SR_delay2(tr - trialRange(1) +1,:)    = sum(SpikesTmp > delay_onset+1 & SpikesTmp < response_onset) ./ (response_onset-(delay_onset+1));
        SR_postcue(tr - trialRange(1) +1,:)   = sum(SpikesTmp > response_onset) ./ postCueDur;
        
        temp = SpikesTmp(SpikesTmp > sample_onset & SpikesTmp < response_onset);
        spikes(1,tr) = {(temp + T.cueOnset(tr)).*1000};
    end
    
    right_trials = (trialTypeInRange==1 & photoStimTrialInRange==0);
    left_trials =  (trialTypeInRange==2 & photoStimTrialInRange==0);

    unit_info.spikes = spikes;
    unit_info.spikes_presample  = SR_presample;
    unit_info.spikes_delay      = SR_delay;
    unit_info.spikes_postcue    = SR_postcue;
    unit_info.delay_activity    = SR_delay-SR_presample;
    unit_info.response_activity = SR_postcue-SR_presample;
    unit_info.trial_type        = trialTypeInRange;
    unit_info.photo_stim        = photoStimTrialInRange;

    unit_info.delay_selectivity    = mean(SR_delay(right_trials)) - mean(SR_delay(left_trials));
    unit_info.response_selectivity = mean(SR_postcue(right_trials)) - mean(SR_postcue(left_trials));

    SA = SR_sample-SR_presample;
    DA = SR_delay-SR_presample;
    RA = SR_postcue-SR_presample;

    norm_DA = (mean(SR_delay(right_trials)) - mean(SR_delay(left_trials))) ./ ...
              (mean(SR_delay(right_trials)) + mean(SR_delay(left_trials)));
    
              (mean(DA(right_trials)) + mean(DA(left_trials)));
    norm_RA = abs(mean(RA(right_trials)) - mean(RA(left_trials))) ./ ...
              (mean(RA(right_trials)) + mean(RA(left_trials)));

    unit_info.norm_delay_selectivity    = norm_DA;
    unit_info.norm_response_selectivity = norm_RA;

    psth_correct_right = mean(PSTH(trialTypeInRange==1 & photoStimTrialInRange==0,:));
    psth_correct_left  = mean(PSTH(trialTypeInRange==2 & photoStimTrialInRange==0,:));

    unit_info.PSTH     = PSTH(trialTypeInRange<3 & photoStimTrialInRange==0,:);
    unit_info.normPSTH = (unit_info.PSTH-nanmin(unit_info.PSTH,[],2)) ./ ... 
                         (nanmax(unit_info.PSTH,[],2));

    % Cell type based on Spike width
    if unit.SpikeWidth > 0.5
        unit_info.type = "RS";
    elseif unit.SpikeWidth < 0.35
        unit_info.type = "FS";
    else 
        unit_info.type = "NA";
    end
    
    % Depth of recording site
    unit_info.depth = unit.Depth;

    % Determine if delay/response selective 
    [p, h, stats] = ranksum(SR_delay(right_trials), SR_delay(left_trials));
    unit_info.selective = h;
    % % % Determine Delay-selective neurons
    if h && (nanmean(SR_delay(right_trials)) > nanmean(SR_delay(left_trials)))
        unit_info.delay_preference = 1;
        DA_preferred = DA(right_trials);
        SR_delay_preferred = SR_delay(right_trials);
        SR_presample_preferred = SR_presample(right_trials);
        SR_delay1_preferred = SR_delay1(right_trials);
        SR_delay2_preferred = SR_delay2(right_trials);
    elseif h && (nanmean(SR_delay(right_trials)) < nanmean(SR_delay(left_trials)))
        unit_info.delay_preference = 2;
        DA_preferred = DA(left_trials);
        SR_delay_preferred = SR_delay(left_trials);
        SR_presample_preferred = SR_presample(left_trials);   
        SR_delay1_preferred = SR_delay1(left_trials);
        SR_delay2_preferred = SR_delay2(left_trials);
    else
        unit_info.delay_preference = 0;    
        DA_preferred = DA(trialTypeInRange<3 & photoStimTrialInRange==0);
        SR_delay_preferred = SR_delay(trialTypeInRange<3 & photoStimTrialInRange==0);
        SR_presample_preferred = SR_presample(trialTypeInRange<3 & photoStimTrialInRange==0);
        SR_delay1_preferred = SR_delay1(trialTypeInRange<3 & photoStimTrialInRange==0);
        SR_delay2_preferred = SR_delay2(trialTypeInRange<3 & photoStimTrialInRange==0);
    end

    % % % Determine Ramping-up or Ramping-down neurons
    [~,h,stats] = signrank(DA_preferred);
    unit_info.delay_sig = h;
    if h && (nanmean(SR_delay_preferred) > nanmean(SR_presample_preferred))
        unit_info.ramp_direction = 1;
    elseif h && (nanmean(SR_delay_preferred) < nanmean(SR_presample_preferred))
        unit_info.ramp_direction = -1;
    else
        unit_info.ramp_direction = 0;
    end
end
