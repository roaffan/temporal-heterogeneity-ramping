
%% Load Data:

load_dir  = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\";
save_dir  = load_dir;

load(fullfile(load_dir, 'ramp_units_mlefit.mat'));
load(fullfile(load_dir, 'alm_dataset5_3150ms.mat'));
load(fullfile(load_dir, 'unit_profile.mat'));

%% Get index of unit types:

r_thresh = 0.5;
rs_ramp  = ramp_RS' & ramp_R2 > 0.5;

R2      = [ramp_R2];
cell_id = [ramp_id'];
temp1   = [ramp_tau];

cell_id = cell_id(R2 > r_thresh);
temp1   = temp1(R2 > r_thresh);

tau = nan(755,1);
tau(cell_id) = temp1;

%%% get opto trials:
data.opto = zeros(size(data.conditions,1),...
                  size(data.conditions,2));
for i = 1:size(data.conditions,1)
    temp = profile(i).photo_stim';
    data.opto(i,1:size(temp,2)) = temp>0;
end
%%%

ramp_units = data;
data = create_specific_class_data(ramp_units,cell_id,tau(cell_id));
save(strcat(save_dir,'upward_neurons.mat'),'data','-v7.3');

clearvars -except data ramp_units load_dir save_dir

%% Save sub-populations of ramp units:

% [~,ind] = sort([data.tau]);
% unit_id = [data.cell_id(ind)];
% tau = [data.tau(ind)];
% N = fix(length(unit_id)./3);
% 
% sort_tau(1).id = [1:N]';
% sort_tau(2).id = [N+1:N*2]';
% sort_tau(3).id = [(N*2)+1:N*3]';
% 
% for i = 1:3
%     clear data
%     data = create_specific_class_data(ramp_units, unit_id([sort_tau(i).id]),tau([sort_tau(i).id]));
%     save(strcat(save_dir,['upward_tau',num2str(i),'_neurons.mat']),'data','-v7.3');
% end

%% Extract specific functional cell types

function newdata = create_specific_class_data(data,cell_id,tau)
    newdata.cell_id          = cell_id;
    newdata.spikes           = data.spikes(cell_id,:);
    newdata.conditions       = data.conditions(cell_id,:);
    newdata.opto             = data.opto(cell_id,:);   
    newdata.number_of_trials = data.number_of_trials(cell_id);
    newdata.trial_length     = data.trial_length;
    newdata.tPreCue          = data.tPreCue;
    newdata.time_axis        = data.time_axis;
    newdata.tau              = tau;
end
