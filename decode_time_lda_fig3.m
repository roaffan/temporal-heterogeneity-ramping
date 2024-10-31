%{
In this analysis, spiking data is binned and sorted then used to train
and test a Linear Discriminate Classifier
%50, 100, 250, 500
bins_all=bins_max-glob_itr; %11, 22, 55 110
variables that can be set to adjust the performance of this
analysis include:
bin_size_constant
    Bin size might be set to 50, 100, 250, 500 to get 110, 55, 22, 11 bins 
    (respectively).  It is currently set to 250 to get reasonable 
    classification accuracy while still having a reasonable number of bins. 
permutations
    Controls how many permutations are used in the permutation testing of
    classifier performance.  1000 permutations provides a strong statistical strength,
    but is slower to run.  
trials
    number of trials used for training classifier and the number of trials
    used for testing the classifier.  A higher number does get slightly
    better performance, but also runs much slower.
Other variable configuration options not used in the final implementation,
these configuration options are not fully bug-tested and may be incomplete: 
c_bins, f_bins
    bins used for classification category (c_bins) and bins used as input 
    (f_bins) can be set seperately although there is little use case for 
    this.  Currently, they are set to be the same value and set
    using bin_size_constant
traingroup_size, testgroup_size
    trials can be clustered together to improve classifier performance,
    This option is currently not used (a cluster of 1)
left_ol, right_overlap
    Sets the training and testing bins to overlap on the left and right
    sides.  This introduces a mathematical problem with the independence of
    the bins and was thus not used (no overlap, i.e. overlap set to 0).
glob_itr
    Not used for this analysis.  Used for supplementary figure S5c.  Used to
    cut earlier time bins.
%}

load_dir  = 'C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data\categorized\';
load_file = 'rampup_neurons.mat';

bin_size_constant = 100;

color_resolution = 6;
cmp = getPyPlot_cMap('BuPu', color_resolution, [], '"C:\Users\roaffan\Anaconda3\envs\caiman\python.exe"');
color = lighten(cmp(end-2,:)+[-0.15 0 0.009],1);
colors(1,:) = lighten(cmp(end-1,:)+[-0.15 0 0.009],1);
colors(2,:) = lighten(cmp(end-2,:)+[-0.15 0 0.009],1);
colors(3,:) = lighten(cmp(end-3,:)+[-0.15 0 0.009],1);

load(strcat(load_dir,load_file))

pre = 0;
glob_itr = 0.000; %0, .500, or 1.150
bins_max = (data.trial_length-((glob_itr*1000)+1))/bin_size_constant;

clearvars -except itr chi2 decoding_error...
           bins_max glob_itr bin_size_constant...
           load_dir load_file pre colors

%number of classifier category bins
%50, 100, 250, 500
bins_all = bins_max - glob_itr; %11, 22, 55 110
c_bins=round(bins_all);

%bins of firing rate
f_bins=round(bins_all);

%%load the data and average firing rate in bins for training trials
%randomly subsample without replacement training trials and testing trials
load(strcat(load_dir,load_file))
[unit_count,trial_count]=size(data.spikes);

it=glob_itr*1000; %initial time cut out in sec
f_bin_size=bin_size_constant;
c_bin_size=bin_size_constant;

trials=125;
traingroup_size=1;
testgroup_size=1;

fr_all_train=zeros(f_bins*trials,unit_count);
fr_all_test=zeros(f_bins*trials,unit_count);

left_ol=1;
right_ol=left_ol-1;

%% Bin the spiking rate data and sort into training and testing trials

for trial=1:trials
    unit_index=0;
    for unit=1:unit_count
        unit_index=unit_index+1;
        for bin=1:f_bins
            train_trial=randsample(1:2:data.number_of_trials(unit),traingroup_size);   %train odd   
            test_trial=randsample(2:2:data.number_of_trials(unit),testgroup_size);    % test even
            spikes_count_train=0;
            spikes_count_test=0;
            for tt=train_trial
            spikes_count_train=spikes_count_train+sum(((bin-left_ol)*f_bin_size+it) < data.spikes{unit,tt} & data.spikes{unit,tt} < ((bin+right_ol)*f_bin_size)+it);
            end
            for tt=test_trial
            spikes_count_test=spikes_count_test+sum(((bin-left_ol)*f_bin_size+it) < data.spikes{unit,tt} & data.spikes{unit,tt} < ((bin+right_ol)*f_bin_size)+it);
            end
            fr_all_train((trial-1)*f_bins+bin,unit_index)=(spikes_count_train)/((min(bin,left_ol)+min(f_bins-bin,right_ol))*f_bin_size)/traingroup_size+rand*10^-13;
            fr_all_test((trial-1)*f_bins+bin,unit_index)=(spikes_count_test)/((min(bin,left_ol)+min(f_bins-bin,right_ol))*f_bin_size)/testgroup_size+rand*10^-13;
        end
    end
end

%set the true category according to the bin count
%works best if f_bins is multiple of c_bins
cat_ind_train(f_bins*trials)=0;
for trial=1:trials
    for bin=1:f_bins
        cat_ind_train((trial-1)*f_bins+bin)=ceil((bin/f_bins)*c_bins);
    end
end

%% run the classifier.  See the matlab help on classify for options

%'linear' uses a linear discriminate analysis
[result,err,posterior,logp,coef] = classify(fr_all_test,fr_all_train,categorical(cat_ind_train)','linear');

%% Analyze and plot the classifier performance
posterior_averaged=zeros(f_bins,c_bins);
expected=ones(f_bins,c_bins)/f_bins;
DoF=f_bins*(c_bins-1);

%Intialize values
rmean(f_bins)=0;
abserrormean(f_bins)=0;
abserrorstd(f_bins)=0;
for bin=1:f_bins
    posterior_averaged(bin,:)=mean(posterior(bin:f_bins:(f_bins*(trials-1)+bin),:),1);
    rerror(bin)=std(grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin))));
    rmean(bin)=mean(grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin))));
    abserrormean(bin)=mean(abs(cat_ind_train(bin)-grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin)))));
    abserrorstd(bin)=std(abs(cat_ind_train(bin)-grp2idx(result(bin:f_bins:(f_bins*(trials-1)+bin)))));
end

color_resolution = 50;
mymap = getPyPlot_cMap('BuPu', color_resolution, [], '"C:\Users\roaffan\Anaconda3\envs\caiman\python.exe"');
mymap = lighten(mymap,.9);

figure;
hold on;
imagesc(log(1/(100)+posterior_averaged),[log(.02) log(.5)]); colormap(mymap);
yl1=ylim;
xl1=xlim;

for bin_i=1:c_bins
   max_f_bin=find(max(posterior_averaged(:,bin_i))==posterior_averaged(:,bin_i));
   plot(bin_i,max_f_bin,'r.','MarkerSize',12,'Color',mymap(50,:));
end
xlim([.5 bins_max])
ylim([.5 bins_max])
xlabel('Actual time (s)')

colorbar

ylabel('Decoded time (s)')
set(gca,'ytick',[pre/c_bin_size,(1150+pre)/c_bin_size,(3150+pre)/c_bin_size])
set(gca,'yticklabel',{[-3.15,-2,0]})

set(gca,'xtick',[pre/c_bin_size,(1150+pre)/c_bin_size,(3150+pre)/c_bin_size])
set(gca,'xticklabel',{[-3.15,-2,0]})
pbaspect([1 1 1])

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

fig_dir = "C:\Users\roaffan\OneDrive\Projects\alm_project_2024\figures_315s\";
fig_file = "lda_mtx_315s";
panel_path = fullfile(fig_dir, fig_file);
print(panel_path, '-dpdf', '-r300')

%% final figure

mdl = fitlm(.001*f_bin_size*(cat_ind_train(round(1150./c_bin_size):round(3150./c_bin_size))-1),...
       .001*c_bin_size*abserrormean(round(1150./c_bin_size):round(3150./c_bin_size)),...
       'y ~ x1');

figure;

hold on 
plot(.001*f_bin_size*(cat_ind_train(round(1150./c_bin_size):round(3150./c_bin_size))-1),...
     mdl.Coefficients{1,1}+mdl.Coefficients{2,1}*.001*f_bin_size*cat_ind_train(round(1150./c_bin_size):round(3150./c_bin_size)),...
     'r-','Color',colors(2,:))
scr=scatter((.001*f_bin_size*(cat_ind_train(round(1150./c_bin_size):round(3150./c_bin_size))-1)),...
        .001*c_bin_size*abserrormean(round(1150./c_bin_size):round(3150./c_bin_size)),12,'filled');
scr.MarkerFaceColor = colors(2,:);

ylim([0 .5]);

xlabel('Time (s)')
ylabel('Error in decoding (s)')

set(gca,'xtick',[0, 1.15, 3.15])
set(gca,'xticklabel',{[-3.15, -2, 0]})
set(gca,'fontname','Arial','fontsize',10)

set(gcf,'color','w')
set(gcf,'paperunits','inches','units','inches')
set(gcf,'pos',[3 2.5 3.11 2.20])
set(gcf,'PaperPositionMode','auto')

axis square

%% Estimate decoding error

diagonal_bins(f_bins)=0;
for i=1:f_bins
    diagonal_bins(i)=posterior_averaged(i,i);
end

overall_averaged_error = mean(.001*c_bin_size*abserrormean);
decoding_error      = .001*c_bin_size*abserrormean(1:31);

%%
disp("saving...")
save_dir  = 'C:\Users\roaffan\OneDrive\Projects\alm_project_2024\data\categorized\';
save_file = ['lda_results.mat'];
save(strcat(save_dir, save_file),...
     'decoding_error','-v7.3');

