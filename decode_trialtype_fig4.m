
%% Load all the data

load_dir   = 'C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data\categorized\';
load_file  = 'rampup_neurons.mat';
total_time = 3150;
         
%% Decode trial type using population activity

clear spikes_trial_ind_all spikes_trial_ind_tmp spikes_trial_ind_mean spikes_trial_ind_std...
    spikes_trial_ind_all_lda_bin_z spikes_trial_ind_all_lda_full_bin_z colors mymap
time_bin = 100;

load(strcat(load_dir,load_file))

n_neurons = length(data.cell_id);
for cell_index=1:n_neurons
%     cell_index = unit_list(cell_index);
    for cat_no = 1:2
        if cat_no==1 % Sensory cue for Right(6kHz)
            subtrial_indexes = find(data.conditions(cell_index,:)==1 | data.conditions(cell_index,:)==3);
        elseif cat_no==2 % Sensory cue for Left(12kHz)
            subtrial_indexes = find(data.conditions(cell_index,:)==2 | data.conditions(cell_index,:)==4);
        end
        clear spikes_trial_ind_tmp
        for trial=1:length(subtrial_indexes)
            subtrial = subtrial_indexes(trial);
            spikes_trial_ind = data.spikes{cell_index,subtrial};

            % truncate decimals from spike index
            spikes_trial_ind = fix(spikes_trial_ind);
            for bin_ind = 1:total_time/time_bin
                spikes_trial_ind_tmp(trial,bin_ind) = ...
                    length(spikes_trial_ind(find(spikes_trial_ind>=time_bin*(bin_ind-1)+1 & spikes_trial_ind<=time_bin*bin_ind)));
                spikes_trial_ind_all{cell_index,cat_no,bin_ind,trial} = spikes_trial_ind_tmp(trial,bin_ind);
            end
        end
        for bin_ind = 1:total_time/time_bin
            spikes_trial_ind_mean{cell_index,cat_no,bin_ind} = mean(spikes_trial_ind_tmp(:,bin_ind));
            spikes_trial_ind_std{cell_index,cat_no,bin_ind} = ...
                std(spikes_trial_ind_tmp(:,bin_ind))/sqrt(length((spikes_trial_ind_tmp(:,bin_ind))));
        end
    end
end

% concatinate trial and category data and convert cell to matrix for LDA
clear spikes_trial_ind_all_lda_full_bin...
    cat_ind_full_bin spikes_trial_ind_all_lda_bin cat_ind_bin

for bin_ind = 1:size(spikes_trial_ind_all,3)
    bin_ind = bin_ind
    spikes_trial_ind_all_lda = [];
    cat_ind = [];
    for cat_no = 1:2
        spikes_trial_ind_all_tmp = squeeze(spikes_trial_ind_all(:,cat_no,bin_ind,:));
        [ii, jj] = find(cellfun(@isempty,spikes_trial_ind_all_tmp));
        spikes_trial_ind_all_tmp = spikes_trial_ind_all_tmp(:,1:min(jj)-1); % eliminate empty trials
        spikes_trial_ind_all_tmp = cell2mat(spikes_trial_ind_all_tmp);
        spikes_trial_ind_all_lda = [spikes_trial_ind_all_lda spikes_trial_ind_all_tmp];
        cat_ind = [cat_ind cat_no*ones(1,size(spikes_trial_ind_all_tmp,2))]; % this is a bit redundant
    end
    spikes_trial_ind_all_lda_bin{bin_ind,:,:} = spikes_trial_ind_all_lda;
end

disp("Performing LDA Classification...")
clear corr_bin_inter incorr_bin_inter same_cat_set_iter diff_cat_set_iter
i = 1;
nn = 20; %number of iterations temp at 2 to run faster to test
n_trial = size(spikes_trial_ind_all_lda_bin{1},2);
tic
while i<= nn
    disp(strcat("iteration #",num2str(i)))
    rand_vec   = randperm(n_trial); % randomly sample trials
    rand_units = randperm(n_neurons-1); % randomly sample cells
    clear corr_bin incorr_bin same_cat_set diff_cat_set
    for bin_ind = 1:size(spikes_trial_ind_all,3)
        bin_ind = bin_ind
        train_len = round(n_trial*0.8);
        test_len = n_trial - train_len;
        cat_ind_train = cat_ind(rand_vec(1:train_len));
        cat_ind_test = cat_ind(rand_vec(train_len+1:end));
        for bin_ind_test = 1:size(spikes_trial_ind_all,3)             
            spikes_trial_ind_all_lda_test = spikes_trial_ind_all_lda_bin{bin_ind_test}(rand_units,rand_vec(train_len+1:end));
            spikes_trial_ind_all_lda_train = spikes_trial_ind_all_lda_bin{bin_ind}(rand_units,rand_vec(1:train_len));
            
            % reduce the dimensionality to full rank 
            [~,colind] = rref(spikes_trial_ind_all_lda_test');
            spikes_trial_ind_all_lda_test = spikes_trial_ind_all_lda_test(colind, :);
            spikes_trial_ind_all_lda_train = spikes_trial_ind_all_lda_train(colind, :);
            
            % reduce the dimensionality to full rank 
            [~,colind] = rref(spikes_trial_ind_all_lda_train');
            spikes_trial_ind_all_lda_test = spikes_trial_ind_all_lda_test(colind, :);
            spikes_trial_ind_all_lda_train = spikes_trial_ind_all_lda_train(colind, :);
            
            %%% classify
            [result,err,posterior,logp,coef] = classify(spikes_trial_ind_all_lda_test',...
                                                        spikes_trial_ind_all_lda_train',...
                                                        categorical(cat_ind_train));  

            cat_diff = grp2idx(result)-cat_ind_test';
            corr_bin(bin_ind,bin_ind_test) = length(find(cat_diff == 0));
            incorr_bin(bin_ind,bin_ind_test) = length(find(cat_diff ~= 0));
            for cat_no = 1:4
                res{cat_no,:} = find(grp2idx(result)==cat_no);
                real_data{cat_no,:} = find(cat_ind_test'==cat_no);
            end 
        end
    end
    corr_bin_inter(i,:,:) = corr_bin;
    incorr_bin_inter(i,:,:) = incorr_bin;
    i = i + 1;
end
toc

%% Save decoding results

tic
disp("saving...")
save_dir  = load_dir;
save_file = ['tau',num2str(itr),'ms_','100ms_20n_LOOV.mat'];
save(strcat(save_dir, save_file),'-v7.3');
toc

%% Set up colors

color_resolution = 6;
cmp = getPyPlot_cMap('BuPu', color_resolution, [], '"C:\Users\roaffan\Anaconda3\envs\caiman\python.exe"');
color = lighten(cmp(end-2,:)+[-0.15 0 0.009],1);
colors(1,:) = lighten(cmp(end-1,:)+[-0.15 0 0.009],1);
colors(2,:) = lighten(cmp(end-2,:)+[-0.15 0 0.009],1);
colors(3,:) = lighten(cmp(end-3,:)+[-0.15 0 0.009],1);
color_resolution = 50;
mymap = getPyPlot_cMap('magma', color_resolution, [], '"C:\Users\roaffan\Anaconda3\envs\caiman\python.exe"');
mymap = lighten(mymap,.9);

%% plot - heatmap for decoding accuracy 

font_size = 10;
lda_mtx = squeeze(mean(corr_bin_inter(:,:,:)))'./...
     (squeeze(mean(corr_bin_inter(:,:,:)))'+squeeze(mean(incorr_bin_inter(:,:,:)))')*100;

figure;imagesc((lda_mtx(1:31,1:31)));
colormap(mymap);colorbar;
caxis([50,100]);
set(gca, 'YDir', 'reverse');

cbh=colorbar();
set(cbh,'Ytick',[50:10:100])
set(get(cbh,'label'),'string','Accuracy (%)');

set(gca,'xtick',[1, 11.5, 31.5])
set(gca,'xticklabel',{'-3.15','-2','0'})
set(gca,'ytick',[1, 11.5, 31.5])
set(gca,'yticklabel',{'-3.15','-2','0'})
xlabel('Train time (s)');
ylabel('Test time (s)');

% % hold on, contour(lda_mtx,[56.9 56.9],'--k','Linewidth',1);
hold on, contour(lda_mtx,[70.4 70.4],'--k','Linewidth',1);
hold on, contour(lda_mtx,[84.8 84.8],'-.k','Linewidth',1);
hold on, contour(lda_mtx,[94.8 94.8],':k','Linewidth',1);

yRange = ylim();
hold on; plot([11.5 11.5] ,yRange, 'w-','LineWidth',1.75)
hold on; plot(yRange, [11.5 11.5], 'w-','LineWidth',1.75)
set(gca,'fontname','Arial','fontsize',10)
axis equal
ylim([11.5, 31.5]), xlim([11.5, 31.5]);
set(gcf,'color','w')
set(gcf,'paperunits','inches','units','inches')
set(gcf,'pos',[3 2.5 2.7014    1.8473])
set(gcf,'PaperPositionMode','auto')

fig_dir = "G:\Shared drives\alm_project\figures_315s\";
fig_file = "ctc_mtx_315s_tau"+num2str(itr);
panel_path = fullfile(fig_dir, fig_file);
print(panel_path, '-dpdf', '-r300')
