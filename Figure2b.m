%%%Get behavior counts and calculate numerical difference pre vs post

clear
clc
close all

%%set up folders
base_path = pwd;
base_excel_dir = 'BorisOutput';
base_contour_dir = 'VocContours';
save_dir = fullfile('OutputPlots\BehavAnalysis');

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

scatter_colors = [1,0,0;0,0.800000000000000,0;
    0,0.400000000000000,0.800000000000000;
    0.800000000000000,0,0.800000000000000;
    1,0.200000000000000,0.600000000000000;
    1,0.501960784313726,0;
    0.627450980392157,0.627450980392157,0.627450980392157];

excel_files = dir(fullfile(base_excel_dir, '*_reorg*'));
excel_fnames = {excel_files.name}';

%%in case behaviors differ across files, let's make a master list
behavs = {'Male follows female', 'Fight', 'Male grooms himself'};

behav_labels = {'MFF', 'Fight', 'MG'};

if ~exist(save_dir)
    mkdir(save_dir);
end

%%get the sinFM data
%%load in data
load('OrganizedVocStruc.mat');
load('VocData-tsne.mat');

%%get pairing info
num_cols = size(InfoStruc,2);
for i = 1:size(InfoStruc,1)
    phase = AllData(strcmp({AllData.AudioFname}, InfoStruc{i,1})).Phase;
    InfoStruc{i,num_cols+1} = phase;
end

prepair_names = {'ch1T0000626', ...
    'ch1T0000627', 'ch1T0000628', 'ch1T0000629', 'ch1T0000642', ...
    'ch1T0000643', 'ch1T0000644'};

postpair_names = {'ch1T0000638', 'ch1T0000639', 'ch1T0000640', ...
    'ch1T0000641', 'ch1T0000677', 'ch1T0000678', 'ch1T0000679'};

var_labels = {'Duration', 'Amp', 'Freq', 'OnsetFreq', 'Slope', 'Phi'};

num_files = length(prepair_names)+length(postpair_names);

%%sort data into pre and postpair (pre and post cohab)
data_prepair = [];
data_postpair = [];
y_data_prepair = [];
y_data_postpair = [];
for i = 1:length(prepair_names)
    idx_tmp = find(strcmp(InfoStruc(:,1), prepair_names{i}) == 1);
    data_prepair = cat(1,data_prepair,InfoStruc{idx_tmp,5});
    y_data_prepair = cat(1,y_data_prepair,InfoStruc{idx_tmp,6});
    cell_data_prepair{i,1} = InfoStruc{idx_tmp,5};
    cell_data_prepair{i,2} = InfoStruc{idx_tmp,6};
end

for i = 1:length(postpair_names)
    idx_tmp = find(strcmp(InfoStruc(:,1), postpair_names{i}) == 1);
    data_postpair = cat(1,data_postpair,InfoStruc{idx_tmp,5});
    y_data_postpair = cat(1,y_data_postpair,InfoStruc{idx_tmp,6});
    cell_data_postpair{i,1} = InfoStruc{idx_tmp,5};
    cell_data_postpair{i,2} = InfoStruc{idx_tmp,6};
end

%%get boundaries of tsne plot
y_data_prepair = cell2mat(cell_data_prepair(:,2));
y_data_postpair = cell2mat(cell_data_postpair(:,2));
all_y_data = cat(1,y_data_prepair,y_data_postpair);
min_val = nanmin(nanmin(all_y_data));
max_val = nanmax(nanmax(all_y_data));

%%I'm just gonna manually put things together instead of reformatting stuff
%%from different sources
date_info = {'0525', '626'; '0525', '627'; '0525', '628'; '0525', '629'; ...
    '0617', '642'; '0617', '643'; '0617', '644'; '0603', '638'; '0603', '639'; ...
    '0603', '640'; '0603', '641'; '0626', '677'; '0626', '678'; '0626', '679'};

%%now, let's try to coordinate these items.
%%I misnamed the video files when I gave them to Jenny, so let's account
%%for that
%%let's get the video dates from the excel files
fdates = {};
all_indices = [];
for i = 1:length(excel_fnames)
    fname_tmp = excel_fnames{i};
    
    %%load in the behavior data
    [~,~,behav_data] = ...
        xlsread(fullfile(base_excel_dir, fname_tmp));
    
    %%now let's coordinate the segments with the behaviors
    for j = 1:length(behavs)
        idx_tmp = find(contains(behav_data(:,7), behavs{j}));
        num_behavs_tmp = length(idx_tmp);
        dur_tmp = [];
        
        for k = 1:length(idx_tmp)
            start_time_tmp = behav_data{idx_tmp(k), 11};
            end_time_tmp = behav_data{idx_tmp(k), 12};
            dur_tmp = cat(1,dur_tmp,end_time_tmp - start_time_tmp);
            
            clear start_time_tmp end_time_tmp voc_raster_tmp voc_times_tmp
        end
        AllDurs{i,j} = dur_tmp;
        AllNumBehavs(i,j) = num_behavs_tmp;
        clear dur_tmp num_behavs_tmp
    end
    clear idx_tmp idx_using audio_fname SegmentData behav_data
end

idx(:,1) = [7:12]; %%prepair
idx(:,2) = [1:6]; %%postpair

%%let's plot the numerical data
figure('color','w'); hold on;
%%sort by pre/post
% AllNumBehavs = cell2mat(AllNumBehavs);
for i = 1:size(AllNumBehavs, 2) %%for each behavior
    avg_behav_num_pre = nanmean(AllNumBehavs(idx(:,1),i));
    avg_behav_num_post = nanmean(AllNumBehavs(idx(:,2),i));
    std_pre = nanstd(AllNumBehavs(idx(:,1),i));
    std_post = nanstd(AllNumBehavs(idx(:,2),i));

    bar([i-.2 i+.2], [avg_behav_num_pre avg_behav_num_post], 'facecolor', 'k');
    line([i-.2 i-.2], [avg_behav_num_pre avg_behav_num_pre+std_pre], 'color', 'k', 'linew', 2);
    line([i-.3 i-.1], [avg_behav_num_pre+std_pre avg_behav_num_pre+std_pre], 'color', 'k', 'linew', 2);
    line([i+.2 i+.2], [avg_behav_num_post avg_behav_num_pre+std_post], 'color', 'k', 'linew', 2);
    line([i+.1 i+.3], [avg_behav_num_pre+std_post avg_behav_num_pre+std_post], 'color', 'k', 'linew', 2);

    for j = 1:length(idx) %%for each rec
        data_tmp_pre = AllNumBehavs(idx(j,1),i);
        data_tmp_post = AllNumBehavs(idx(j,2),i);
        line([i-.2 i+.2], [data_tmp_pre data_tmp_post], 'color', [.6 .6 .6], 'linew', 1);
        scatter(i-.2, data_tmp_pre, 40, 'filled', 'markerfacecolor', scatter_colors(j,:));
        scatter(i+.2, data_tmp_post, 40, 'filled', 'markerfacecolor', scatter_colors(j,:));
    end
end
set(gca, 'xtick', ([1:1:size(AllNumBehavs,2)]));
set(gca, 'xticklabel', behav_labels);
xlabel('Behavior');
ylabel('Count');

%%let's run the paired ttests
for i = 1:size(AllNumBehavs,2) %%for each behavior
    [h,p(i,1),stats{i}] = ttest2(AllNumBehavs(idx(:,1),i), AllNumBehavs(idx(:,2),i));
end

text(1, 75, sprintf('p=%.2f', p(1)));
text(1, 70, sprintf('p=%.2f', p(2)));
text(1, 65, sprintf('p=%.2f', p(3)));

saveas(gcf, fullfile(save_dir, 'BehavCounts.jpg'));
saveas(gcf, fullfile(save_dir, 'BehavCounts.svg'));
