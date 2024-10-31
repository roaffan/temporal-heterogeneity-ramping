function fit_params_exG(ii)

%Loads up all the initialization stuff
initiate_ML_exG

%Checks if job was submitted to cluster
if params.submit 
    %Turns string of cell number into double
    if isstring(ii) || ischar(ii)
        cell_no=str2double(ii);
    else
        cell_no=ii;
    end
else
    %Pulls cell number
    cell_no=ii;
end

% % %Take care of missing files
% unit_id = [6,22,110,184];
% ii = unit_id(ii);
% cell_no = ii;

%Loads the spike data
load(params.data);
load(params.index);

%Exit if not pyramidal or fast-spiking:
if ~index.pyramidal(ii)
    return
end

%Exit if not significantly ramp-up/down:
if (index.ramp_down(ii) | index.ramp_up(ii))==0
    return
end

%Variables that need to be passable between functions
global f_spikes f_spikes_evenodd t v c train_test_flag with_T lb ub x trial_length st ut y

%particle swarming options
hybridopts = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000);
options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',{@fmincon,hybridopts}); 

%Gets length of trial
trial_length = data.trial_length;
%Gets time before trial
time_before = data.time_before+1150;
%Gets time after go cue
time_after = data.time_after;

trial_length = trial_length-(time_before+time_after);

%Holds lower bounds
lb = [];
%Holds upper bounds
ub = [];
%Prints cell number
cell_no = cell_no
%Holds all of the spikes
f_spikes = [];
f_cond   = [];
%loops through the trials to reformat for solver
for j=1:data.number_of_trials(cell_no)
    %Will hold all the spikes in the time window I care about
    %Each millisecond bin either has a spike (1) or doesn't (0)
    f_trial = zeros(trial_length,1);
    %Gets the spikes in the given trial
    temp1 = fix(data.spikes{cell_no,j})-time_before;
    temp1 = temp1(0<temp1 & temp1<trial_length);
    %Only loads up the spikes during the trial
    f_trial(temp1) = 1;

    %Stores the spikes
    f_spikes=[f_spikes; f_trial];
 
    %Will hold all the behavioral/trial condition when spike occured
    f_cond_trl  = zeros(trial_length,1);
    %Gets the spikes in the given trial
    temp2 = data.conditions(cell_no,j);
    f_cond_trl= repmat(temp2,trial_length,1);
    %Stores the trial condition for spikes
    f_cond=[f_cond; f_cond_trl];
end

%Reformats time to repeat the same period over and over again for solver
t = repmat([1:trial_length], 1, length(f_spikes)/trial_length )';
t = t/1000;
t = t-mean(t);

%% Create weights for behavioral conditions:

% correct right
cr=f_cond;
cr(cr==1)=1;
cr(cr==2)=0;

% correct left
cl=f_cond;
cl(cl==1)=0;
cl(cl==2)=1;

v={cr,cl};

%% GLM    

%Loads what equations are going to be used.
fun = @glm_model_exG;

%% Check unit activity type:

if (index.ramp_up(ii) && ~index.trial_type(ii))
    model = 1;
elseif (index.ramp_up(ii) && index.trial_type(ii))
    model = 2;
elseif (index.ramp_down(ii) && ~index.trial_type(ii))
    model = 3;
elseif (index.ramp_down(ii) && index.trial_type(ii))
    model = 4;    
end    

%% Exponential-modified Gaussian:
tic
if model==1 || model==3
    with_T = 1;
elseif model==2 || model==4
    with_T = 2;
end

%number of terms in the model (baseline and EMG)
n = 2;

%Tau lower bound (model uses 1/tau, so 1/0.10 = 10 s)
%Tau upper bound (model uses 1/tau, so 1/100 = 10 ms)
lb.tau = 0.10;
ub.tau = 10;

% %Mu for ramps will be from middle of time window to the go cue
% lb.mu = min(t);
% ub.mu = min(t)+.5;

%Sigma will be very small
lb.sig = 0.001;
ub.sig = 0.1;

if model==1 
    %Baseline 
    lb.o = 10^-5;
    ub.o = 1/n;
    %Amplitude 
    lb.peak = 10^-5;
    ub.peak = 1/n;
    %Stores ramp parameter bounds
    lbR = lb;
    ubR = ub;
    %Runs the solver and stores all the information from the fit.
    [xR, xR_tmp, LL_R, LL_R_tmp, LL_R_trial_2blocks, aic_R, bic_R, R] = run_solver(params, fun, options);
elseif model==2
    %Baseline 
    lb.o = 10^-5;
    ub.o = 1/n;
    %Amplitude 
    lb.peak = 10^-5;
    ub.peak = 1/n;
    lb.peak2 = 10^-5;
    ub.peak2 = 1/n;
    % Trial-type Specific Ramp-up
    lbRS = lb;
    ubRS = ub;
    [xRS, xRS_tmp, LL_RS, LL_RS_tmp, LL_RS_trial_2blocks, aic_RS, bic_RS, RS] = run_solver(params, fun, options);
elseif model==3
    %Baseline 
    lb.o = xC;
    ub.o = 1/n;
    %Amplitude 
    lb.peak = -xC;
    ub.peak = -(10^-5);
    % Ramp-down/decay
    lbD = lb;
    ubD = ub;
    [xD, xD_tmp, LL_D, LL_D_tmp, LL_D_trial_2blocks, aic_D, bic_D, D] = run_solver(params, fun, options);
elseif model==4
    %Baseline 
    lb.o = xC;
    ub.o = 1/n;
    %Amplitude 
    lb.peak = -xC;
    ub.peak = -(10^-5);
    lb.peak2 = -xC;
    ub.peak2 = -(10^-5);
    % Trial-type Specific Decay
    lbDS = lb;
    ubDS = ub;
    [xDS, xDS_tmp, LL_DS, LL_DS_tmp, LL_DS_trial_2blocks, aic_DS, bic_DS, DS] = run_solver(params, fun, options);    
end

toc

%%
%Clears out data to make the size of saved file smaller
clear data;
%Saves results of script
save(sprintf('%s/glm_cell_%i_exG.mat', params.folder_results,cell_no));

%Checks if it was submitted on cluster
if params.submit
    %Closes program on cluster
    exit
end

end


function [xC, xC_tmp, LL_C, LL_C_tmp, LL_C_trial_2blocks, aic_C, bic_C, C] =...
         run_solver(params, fun, options)

global train_test_flag f_spikes lb ub
%global f_spikes t train_test_flag with_T lb ub trial_length
%create vectors from upper and lower bounds structures since
%particleswarm requires vector input
lb_cell=struct2cell(lb);
lb_vec=[lb_cell{:}];
ub_cell=struct2cell(ub);
ub_vec=[ub_cell{:}];
train_test_flag=1; %training - return the fit sum
tic
stop_loop=0;
c=1;
LL_C_tmp_max=-Inf;
i=0;
while stop_loop==0
   %Helps to handle occasions when probability is not between 0 and 1.
   try
       i=i+1;
       xC_tmp(i,:)=particleswarm(fun,length(lb_vec),lb_vec,ub_vec,options);
       LL_C_tmp{i}=-fun(xC_tmp(i,:));
       if LL_C_tmp_max>=LL_C_tmp{i}
           c=c+1
       else
           LL_C_tmp_max=LL_C_tmp{i};
           xC_tmp_max=xC_tmp(i,:);
           c=1
       end
       if c>=params.N
           stop_loop=1;
       end
   catch
       c
   end
       
end
LL_C=LL_C_tmp_max;
xC=xC_tmp_max;
if isfield(lb,'t') || length(lb_vec) == 1
    train_test_flag=2; %compute LL per trial - divide trials into 2 blocks
    LL_C_trial_2blocks=-fun(xC);
else
    LL_C_trial_2blocks = [];
end
[aic_C,bic_C]=aicbic(LL_C,length(lb_vec),length(f_spikes));          
toc
train_test_flag=0; %testing - rerun the fit
C=fun(xC);
clear f_spikes t train_test_flag with_T lb ub trial_length
end

%Make sure to compile with matlab/2016a

%compile on the cluster
%mcc -mv -o stim_spec_time_cells_16a ML_fit_params.m -R -singleCompThread -R -nodisplay
