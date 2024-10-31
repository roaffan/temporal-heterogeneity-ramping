%% initialize
%Before running make sure the following is good to go
%Is it set to submit
%What is N at

%What folders need to be added
addpath(genpath('/projectnb/ecog-eeg/alm_cells/'));

mainDir = '/projectnb/ecog-eeg/alm_cells/';

%Where are the results being saved
params.folder_results = strcat(mainDir,'results/exG_3150ms_august23');
params.folder_mat = strcat(mainDir,'data');
%The actual data
params.data = strcat(params.folder_mat,'/alm_data.mat'); 
params.index = strcat(params.folder_mat,'/alm_index.mat'); 

%This is how many times I want to try to best
%the previous fit.
params.N = 20; %20 for final
%params.submit must = 1 if submited as job
params.submit =1;

