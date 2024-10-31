
clear

%% Load data

main_dir = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\";

% folder containing output files of model fitting (one for each neuron)
mle_dir = strcat(main_dir,'data\','3150ms_fit_mle\v6_315s_8-16-23\');

load(fullfile(main_dir,'data','unit_profile.mat'));
unit_data = profile;

load(fullfile(main_dir,'data','alm_index.mat'));
unit_type = index;

%% Iterate through the MLE result files

missing_files = [];
for unit_i = 1:755

    filename = fullfile(mle_dir, strcat('glm_cell_', num2str(unit_i),...
                  '_exG_3150ms.mat'));

    if exist(filename, 'file') ~= 2
        unit_data(unit_i).tau = nan;
        unit_data(unit_i).sigma = nan;
        unit_data(unit_i).mu = nan;
        unit_data(unit_i).fmodel  = [];
        unit_data(unit_i).fspikes = [];    
        
        if unit_type.ramp_up(unit_i) || unit_type.ramp_down(unit_i)
            if unit_type.pyramidal(unit_i) || unit_type.fastspike(unit_i)
                missing_files = [missing_files, unit_i];
            end
        end
    elseif exist(filename, 'file') == 2 
        load(filename);
        % Determine parameter for specific model
        if exist('xR'), x = xR; fmod = R;
        elseif exist('xRS'), x = xRS; fmod = RS;
        elseif exist('xD'),  x = xD; fmod = D;
        elseif exist('xDS'), x = xDS; fmod = DS;
        end
        % Store parameters and results 
        if exist('xR') || exist('xRS')
            unit_data(unit_i).tau = 1./x(1);
            unit_data(unit_i).mu = -x(2) + (3.15/2);
            unit_data(unit_i).sigma = x(3);
        elseif exist('xD') || exist('xDS')
            unit_data(unit_i).tau = 1./x(2);
            unit_data(unit_i).mu = -x(4) + (3.15/2);
            unit_data(unit_i).sigma = x(3);
        end
        unit_data(unit_i).fmodel  = nanmean(reshape(fmod, [trial_length, numel(fmod)/trial_length])',1);
        unit_data(unit_i).fspikes = nanmean(reshape(f_spikes, [trial_length, numel(f_spikes)/trial_length])',1);
    end
    
clearvars -except mle_dir unit_data unit_type missing_files go_cue_sec...
           time_length time_info
end

%% Save data
filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";
save(fullfile(filename, 'alm_mlfits.mat'),'unit_data',...
                        'unit_type','missing_files','-v7.3');

%% Get index for all unit types

unit_id = 1:1:755;

RS_ramp = unit_type.pyramidal & unit_type.ramp_up & [unit_data.delay_sig];
FS_ramp = unit_type.fastspike & unit_type.ramp_up;
RS_decay = unit_type.pyramidal & unit_type.ramp_down & [unit_data.delay_sig];
FS_decay = unit_type.fastspike & unit_type.ramp_down;

time_length = 3150 + 1;

%% Create a matrix of normalized PSTHs

go_cue = find(time_info.tAxisPSTH==0);
delay_onset = go_cue-2000-1150;

for unit_i = 1:755
    temp = nanmean(unit_data(unit_i).PSTH,1);
    PSTHs(unit_i,:) = temp(delay_onset:go_cue);

    toi = unit_data(unit_i).trial_type<3 & unit_data(unit_i).photo_stim==0;
    trial_ind = unit_data(unit_i).trial_type(toi);
    temp = nanmean(unit_data(unit_i).PSTH(trial_ind==1,:),1);
    PSTH_right(unit_i,:) = temp(delay_onset:go_cue);
    temp = nanmean(unit_data(unit_i).PSTH(trial_ind==2,:),1);
    PSTH_left(unit_i,:) = temp(delay_onset:go_cue);

    if isempty(unit_data(unit_i).fmodel)
        fModel(unit_i,:) = nan(1,time_length);
        fSpikes(unit_i,:) = nan(1,time_length);
        normModel(unit_i,:) = nan(1,time_length);
        normPSTH(unit_i,:) = nan(1,time_length);
    else
        temp = unit_data(unit_i).fmodel;
        fModel(unit_i,:) = temp;
        temp = (temp - min(temp));
        temp = temp ./ max(temp);
        normModel(unit_i,:) = temp;   
        
        temp = smooth(unit_data(unit_i).fspikes,101);
        fSpikes(unit_i,:) = temp;
        temp = PSTHs(unit_i,:);
        temp = (temp - min(temp));
        temp = temp ./ max(temp);
        normPSTH(unit_i,:) = temp;   
    end
end

mu    = [unit_data.mu]';
tau   = [unit_data.tau]';
sigma = [unit_data.sigma]';
depth = [unit_data.depth]';
selectivity = [unit_data.norm_delay_selectivity];

%% Estimate goodness of fit between PSTHs and model fit

for unit_i = 1:755
    if isempty(unit_data(unit_i).fspikes)
        Rsquared(unit_i,:) = nan;
    else
        mdl = fitlm(fSpikes(unit_i,:), fModel(unit_i,:));
        Rsquared(unit_i,:) = mdl.Rsquared.Ordinary;
        close;
    end
end

%% Save updated unit data

filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";

save(fullfile(filename, 'alm_mlfit.mat'),'unit_id','unit_type','time_info',...
    'PSTHs','PSTH_right','PSTH_left','normPSTH','normModel',...
    'fSpikes','fModel','mu','tau','sigma','selectivity','Rsquared','-v7.3');

%% Save ramp vs decay unit data

%%% save ramps:
ramp_id = unit_id(RS_ramp);
ramp_RS = RS_ramp(RS_ramp);
ramp_FS = FS_ramp(unit_type.ramp_up);
ramp_fSpikes = fSpikes(RS_ramp,:);
ramp_fModel = fModel(RS_ramp,:);
ramp_normPSTH = normPSTH(RS_ramp,:);
ramp_normModel = normModel(RS_ramp,:);
ramp_mu = mu(RS_ramp);
ramp_tau = tau(RS_ramp);
ramp_sigma = sigma(RS_ramp);
ramp_depth = depth(RS_ramp);
ramp_selectivity = selectivity(RS_ramp);
ramp_R2 = Rsquared(RS_ramp);

filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";
save(fullfile(filename, 'ramp_units_mlefit.mat'),'ramp_id','ramp_RS',...
    'ramp_normPSTH','ramp_normModel','ramp_mu','ramp_tau','ramp_sigma',...
    'ramp_depth','ramp_selectivity','ramp_R2','-v7.3');

%%% save decays:
decay_id = unit_id(RS_decay);
decay_RS = RS_decay(RS_decay);
decay_FS = FS_decay(RS_decay);
decay_normPSTH = normPSTH(RS_decay,:);
decay_normModel = normModel(RS_decay,:);
decay_mu = mu(RS_decay);
decay_depth = depth(RS_decay);
decay_tau = tau(RS_decay);
decay_sigma = sigma(RS_decay);
decay_selectivity = selectivity(RS_decay);
decay_R2 = Rsquared(RS_decay);

filename = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data";
save(fullfile(filename, 'decay_units_mlefit.mat'),'decay_id','decay_RS',...
    'decay_normPSTH','decay_normModel','decay_mu','decay_tau','decay_sigma', ...
    'decay_depth','decay_selectivity','decay_R2','-v7.3');

