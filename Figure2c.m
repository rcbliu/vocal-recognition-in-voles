%%%script look at behavior associated vocalizations by peak in 2A

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

%%organize data by prepair and postpair (precohab and postcohab)
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
    dash = strfind(fname_tmp, '-');
    underscore = strfind(fname_tmp, '_');
    fdates{i,1} = excel_fnames{i}(dash(1)+1:dash(2)-1);
    fnum_tmp = fname_tmp(underscore(1)-1);
    %%separate out pre and post
    %%they're already in chronological order otherwise
    if strcmp(fname_tmp(1:dash(1)-1), 'Postpair')
        grp = 2;
    else
        grp = 1;
    end
    
    idx_tmp = find(contains(date_info(:,1), fdates{i,1}));
    idx_using = idx_tmp(str2num(fnum_tmp));
    all_indices(end+1,1) = idx_using;
    
    %%load in the correct contours
    audio_fname = date_info{idx_using,2};
    try
        load(fullfile(base_contour_dir, ...
            sprintf('SegmentData_ch1T0000%i.mat', str2num(audio_fname))));
    catch
        load(fullfile(base_contour_dir, ...
            sprintf('SegmentData_T0000%i.mat', str2num(audio_fname))));
    end
    
    %%get the relevant tsne data
    audio_fname_tmp = sprintf('ch1T0000%i', str2num(audio_fname));
    try
        fm_idx = find(strcmp(prepair_names, audio_fname_tmp));
        sinfm_data = cell_data_prepair{fm_idx,2};
    catch
        fm_idx = find(strcmp(postpair_names, audio_fname_tmp));
        sinfm_data = cell_data_postpair{fm_idx,2};
    end
    %%load in the behavior data
    [~,~,behav_data] = ...
        xlsread(fullfile(base_excel_dir, fname_tmp));
    
    %%now let's coordinate the segments with the behaviors
    for j = 1:length(behavs)
        idx_tmp = find(contains(behav_data(:,7), behavs{j}));
        num_vocs_tmp = [];
        voc_rate_tmp = [];
        voc_raster_tmp = [];
        voc_data_tmp = {};
        voc_locs_tmp = [];
        dur_tmp = [];
        
        for k = 1:length(idx_tmp)
            start_time_tmp = behav_data{idx_tmp(k), 11};
            end_time_tmp = behav_data{idx_tmp(k), 12};
            
            voc_idx_tmp = find(SegmentData.StartTime >= start_time_tmp & ...
                SegmentData.StartTime < end_time_tmp);
            num_vocs_tmp(k,1) = length(voc_idx_tmp);
            voc_rate_tmp(k,1) = num_vocs_tmp(k,1)/(end_time_tmp - start_time_tmp);
            try
                voc_locs_tmp{k,1} = ...
                    sinfm_data(voc_idx_tmp,:);
            catch
                idx_usable = find(voc_idx_tmp < length(sinfm_data));
                voc_locs_tmp{k,1} = ...
                    sinfm_data(idx_usable,:);
            end
            
            %%let's also save the time info for the vocs so we can use it
            %%later to make a raster plot
            voc_raster_tmp(1,1) = start_time_tmp;
            voc_raster_tmp(1,2) = end_time_tmp;
            if ~isempty(voc_idx_tmp)
                voc_raster_tmp = cat(2,voc_raster_tmp, ...
                    SegmentData.StartTime(voc_idx_tmp)');
                voc_data_tmp{end+1,1} = voc_raster_tmp;
            end
            dur_tmp = cat(1,dur_tmp,end_time_tmp - start_time_tmp);
            
            clear start_time_tmp end_time_tmp voc_raster_tmp voc_times_tmp
        end
        AllNumVocs{i,j} = num_vocs_tmp;
        AllRateVocs{i,j} = voc_rate_tmp;
        BehavVocs{i,j} = voc_data_tmp;
        tSNELocs{i,j} = voc_locs_tmp;
        AllDurs{i,j} = dur_tmp;
        if length(voc_locs_tmp) ~= length(dur_tmp)
            disp('Issue')
        end
        clear num_vocs_tmp idx_tmp voc_rate_tmp voc_data_tmp voc_locs_tmp dur_tmp
    end
    clear idx_tmp idx_using audio_fname SegmentData behav_data
end

%%since we already have the tsne locs, let's label each sound with which
%%peak it comes from
% [edges1,edges2,edges3] = ...
%     fn_label_peak_locs(all_y_data); %%saved output of this (because it
%     takes forever to run) for next section
load('tsne_peak_vertices');
edges1 = tsne_peak_vertices{1};
edges2 = tsne_peak_vertices{2};
edges3 = tsne_peak_vertices{3};

for i = 1:size(tSNELocs,1) %%for each rec
    for j = 1:size(tSNELocs,2) %%for each behavior
        data_using = tSNELocs{i,j};
        for k = 1:length(data_using)
            data_tmp = data_using{k};
            data_tmp(:,3) = zeros(size(data_tmp,1),1);
            [in1,on1] = inpolygon(data_tmp(:,1), data_tmp(:,2),...
                edges1(:,1),edges1(:,2));
            idx_1 = find(in1 == 1);
            data_tmp(idx_1,3) = 1;

            [in2,on2] = inpolygon(data_tmp(:,1), data_tmp(:,2),...
                edges2(:,1),edges2(:,2));
            idx_2 = find(in2 == 1);
            data_tmp(idx_2,3) = 2;

            [in3,on3] = inpolygon(data_tmp(:,1), data_tmp(:,2),...
                edges3(:,1),edges3(:,2));
            idx_3 = find(in3 == 1);
            data_tmp(idx_3,3) = 3;

            data_using{k} = data_tmp;
        end
        tSNELocs{i,j} = data_using;
        clear data_tmp idx_* in* on* data_using
    end
end

%%now let's make plots for voc rates just within the 3 main clusters
for i = 1:size(tSNELocs,2) %%for each behavior
    for j = 1:size(tSNELocs,1) %%for each individual
        %%let's get the vocal rate per peak and per behavior
        voc_data_using = tSNELocs{j,i};
        durations_tmp = AllDurs{j,i};

        if isempty(voc_data_using)
            continue
        end

        for k = 1:length(voc_data_using)
            voc_data_tmp = voc_data_using{k};
            voc_rate1_tmp(k) = length(find(voc_data_tmp(:,3) == 1))/durations_tmp(k);
            voc_rate2_tmp(k) = length(find(voc_data_tmp(:,3) == 2))/durations_tmp(k);
            voc_rate3_tmp(k) = length(find(voc_data_tmp(:,3) == 3))/durations_tmp(k);
        end
        voc_rate1{i,j} = voc_rate1_tmp;
        voc_rate2{i,j} = voc_rate2_tmp;
        voc_rate3{i,j} = voc_rate3_tmp;
    end
end

%%calculate the average vocal rate within each behavior.
%%rows = individuals, columns = behaviors
%%need to reorganize so the corresponding data is together.
%%hard coding
idx(:,1) = [7:12]; %%prepair data
idx(:,2) = [1:6]; %%postpair data

for i = 1:size(AllRateVocs,2) %%for each behavior
    for j = 1:size(AllRateVocs,1)/2 %%for each individual
        data_tmp1 = voc_rate1{i,idx(j,1)};
        data_tmp1 = cat(2,data_tmp1,voc_rate1{i,idx(j,2)});
        if ~isempty(data_tmp1)
            avg_voc_rate(j,i,1) = nanmean(data_tmp1);
        end
        data_tmp2 = voc_rate2{i,idx(j,1)};
        data_tmp2 = cat(2,data_tmp2,voc_rate2{i,idx(j,2)});
        if ~isempty(data_tmp2)
            avg_voc_rate(j,i,2) = nanmean(data_tmp2);
        end
        data_tmp3 = voc_rate3{i,idx(j,1)};
        data_tmp3 = cat(2,data_tmp3,voc_rate3{i,idx(j,2)});
        if ~isempty(data_tmp3)
            avg_voc_rate(j,i,3) = nanmean(data_tmp3);
        end
        clear data_tmp1 data_tmp2 data_tmp3
    end
end

%%Let's start by just doing 3 columns per behavior, peaks1-3, and then
%%separate plots for pre/post
fig1 = figure('color','w'); hold on
for i = 1:size(avg_voc_rate,2)
    figure(fig1); %%do pre data
    locs = [-0.2 0 0.2];
    h = bar([i+locs(1) i+locs(2), i+locs(3)], ...
        [nanmean(avg_voc_rate(:,i,1)) nanmean(avg_voc_rate(:,i,2)) nanmean(avg_voc_rate(:,i,3))], ...
        'barw', .8, 'facecolor', 'k');
    for j = 1:size(avg_voc_rate,1)
        line([i+locs(1) i+locs(2) i+locs(3)], ...
            [avg_voc_rate(j,i,1) avg_voc_rate(j,i,2) avg_voc_rate(j,i,3)], ...
            'color', [.6 .6 .6], 'linew', 1)
        scatter(i+locs(1), avg_voc_rate(j,i,1), 40, 'filled', 'markerfacecolor', scatter_colors(j,:));
        scatter(i+locs(2), avg_voc_rate(j,i,2), 40, 'filled', 'markerfacecolor', scatter_colors(j,:));
        scatter(i+locs(3), avg_voc_rate(j,i,3), 40, 'filled', 'markerfacecolor', scatter_colors(j,:));
    end
end
figure(fig1);
set(gca, 'xtick', [1:1:length(behavs)]);
set(gca, 'xticklabel', behav_labels);
ylabel('Vocal Rate');
title('Vocal Rates Behavs tSNE Peaks');
ylim([0 0.6]);

saveas(gcf, fullfile(save_dir, 'VocsInBehavs-tSNE-CollapseExper.jpg'))
saveas(gcf, fullfile(save_dir, 'VocsInBehavs-tSNE-CollapseExper.svg'))


%%let's reorganize the data for an anova
anova_data = [];
anova_labels = [];
for i = 1:size(avg_voc_rate,2) %%for each behav
    for j = 1:size(avg_voc_rate,1) %%for each individual
        anova_data(end+1,1) = avg_voc_rate(j,i,1); %%Cluster 1
        anova_data(end+1,1) = avg_voc_rate(j,i,2); %%Cluster 2
        anova_data(end+1,1) = avg_voc_rate(j,i,3); %%Cluster 3

        anova_labels(end+1,1) = 1; %%voc cluster
        anova_labels(end,2) = j;   %%id
        anova_labels(end,3) = i;   %%behavior

        anova_labels(end+1,1) = 2; %%voc cluster
        anova_labels(end,2) = j;   %%id
        anova_labels(end,3) = i;   %%behavior

        anova_labels(end+1,1) = 3; %%voc cluster
        anova_labels(end,2) = j;   %%id
        anova_labels(end,3) = i;   %%behavior
    end
end
cluster_labels = anova_labels(:,1);
animal_id_labels = anova_labels(:,2);
behavior_labels = anova_labels(:,3);


[p,tbl,stats] = anovan(anova_data, {cluster_labels, behavior_labels, animal_id_labels}, ...
    'model', 'interaction', 'varnames', {'Cluster', 'Behav',  'ID'});
saveas(gcf, fullfile(save_dir, 'VocsInBehavs-tSNE-CollapseExper-ANOVA.jpg'))
figure('color','w'); 
[c,m,h,gnames] = multcompare(stats)