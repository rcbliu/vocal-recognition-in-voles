%%%characterize the raw acoustic features of adult prairie vole USVs

clear
clc
close all

use_patch = 0; %%plot with a box and whisker plot (1) or just a line for mean +/- s.d.
% save_dir = 'S:\Megan\sinFM_GPU_Output\Output_newTsne\AllData';
save_dir = 'OutputPlots\';
base_path = pwd;
[~,~,data] = xlsread('RecordingList');
%%This data is organized such that each row corresponds to three separate
%%audio files (pretest, prepair, posttest). Read in all of that audio data.

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

scatter_colors = [1,0,0;...
    0,0.800000000000000,0;...
    0,0.400000000000000,0.800000000000000;...
    0.800000000000000,0,0.800000000000000;...
    1,0.200000000000000,0.600000000000000;...
    1,0.501960784313726,0;...
    0.627450980392157,0.627450980392157,0.627450980392157];
    
if ~exist(save_dir)
    mkdir(save_dir);
end

% fields_of_interest = [7 9:12 14 18 21];
% features_of_interest = [3 5 6 7 8 10 12 13 14];
features_of_interest = [3 5:8 11 14 17];

for dataset = 2:size(data,1) %%ignore header row
    tic
    
    for j = 2:3
        if j == 2
            date_tmp = data{dataset,8};
            label_tmp = '30MinMFDyads-prePairing';
            audio_fname = data{dataset,9};
        elseif j == 3
            date_tmp = data{dataset,10};
            label_tmp = '30MinMFDyads-postPairing';
            audio_fname = data{dataset,11};
        end
        backslash = strfind(date_tmp, '/');
        month = date_tmp(1:backslash(1)-1);
        day = date_tmp(backslash(1)+1:backslash(2)-1);
        year = date_tmp(backslash(2)+1:end);
        date_str = sprintf('%02d%02d%04d', str2double(month),str2double(day),str2double(year));
        base_path_tmp = fullfile(base_path, label_tmp, date_str, 'Audio');
        audio_path = base_path_tmp;
        contour_path = fullfile('VocContours\RawAcousticFeats');
        
        if strcmp(date_str, '05232021')%%I messed up this file. Ignore it for now.
            continue
        end
        
        fprintf('Starting %s\n', audio_fname);
        tic
        
        s = load(fullfile(contour_path, sprintf('VocData_%s', audio_fname)));
        
        toc
        VocData = s.VocalizationData; clear s
        
        fields = fieldnames(VocData);
        
        for k = 1:length(features_of_interest)
            avg_vals(dataset-1,k,j-1) = nanmean(VocData.(sprintf('%s', fields{features_of_interest(k)})));
        end
        
        num_usvs(dataset-1,j-1) = length(VocData.StartTime);
        
    end
end

%%now run stats
for i = 1:size(avg_vals,2)
    [h,p,ci,stats] = ttest(avg_vals(:,i,1),avg_vals(:,i,2));
    all_p(i,1) = p;
    all_stats(i,1) = stats.tstat;
    clear h p stats
end

figure('color','w'); hold on
scatter_colors = [255 0 0; 0 204 0; 0 102 204; 204 0 204; 255 51 153; 255 128 0; 160 160 160]./255;
for i = 1:size(avg_vals,2)
    for j = 1:size(avg_vals,1)
        line([i-.25 i+.25], [avg_vals(j,i,1) avg_vals(j,i,2)], 'color','k');
        scatter([i-.25 i+.25]', [avg_vals(j,i,1) avg_vals(j,i,2)]', 30, repmat(scatter_colors(j,:),2,1), 'filled');
    end
end
set(gca, 'yscale', 'log');

set(gca, 'xtick', [1:1:size(avg_vals,2)], 'xticklabel', fields(features_of_interest));
saveas(gcf, fullfile(save_dir, 'RawAcousticFeats-Vocs.jpg'));
saveas(gcf, fullfile(save_dir, 'RawAcousticFeats-Vocs.svg'));


%%Also make individual plots for each feature
for i = 1:size(avg_vals,2)
    figure('color','w','position',[360 406 168 212]); hold on
    mean_val_pre = nanmean(avg_vals(:,i,1));
    sd_val_pre = nanstd(avg_vals(:,i,1))./size(avg_vals,1);
    mean_val_post = nanmean(avg_vals(:,i,2));
    sd_val_post = nanstd(avg_vals(:,i,2))./size(avg_vals,1);
    if use_patch == 1
        patch([.7 .7 1.3 1.3], [mean_val_pre+sd_val_pre mean_val_pre-sd_val_pre mean_val_pre-sd_val_pre mean_val_pre+sd_val_pre], ...
            [253,245,230]./255);
        patch([1.7 1.7 2.3 2.3], [mean_val_post+sd_val_post mean_val_post-sd_val_post mean_val_post-sd_val_post mean_val_post+sd_val_post], ...
            [253,245,230]./255);
        line([.7 1.3], [mean_val_pre mean_val_pre], 'color','k','linew',3);
        line([1.7 2.3], [mean_val_post mean_val_post], 'color','k','linew',3);
    end
    for j = 1:size(avg_vals,1)
        line([1 2], [avg_vals(j,i,1) avg_vals(j,i,2)], 'color',[150 150 150]./255, 'linew', 1);
%         scatter([1 2]', [avg_vals(j,i,1) avg_vals(j,i,2)]', 60, repmat(scatter_colors(j,:),2,1), 'filled');
        scatter([1 2]', [avg_vals(j,i,1) avg_vals(j,i,2)]', 30, repmat(scatter_colors(j,:),2,1), 'filled');
    end
    if use_patch == 0
        line([0.9 1.1], [mean_val_pre mean_val_pre], 'color', 'k', 'linew', 2);
        line([1 1], [mean_val_pre+sd_val_pre mean_val_pre-sd_val_pre], 'color','k','linew',1);
        line([.95 1.05], [mean_val_pre+sd_val_pre mean_val_pre+sd_val_pre], 'color','k','linew',1);
        line([.95 1.05], [mean_val_pre-sd_val_pre mean_val_pre-sd_val_pre], 'color','k','linew',1);
        
        line([1.9 2.1], [mean_val_post mean_val_post], 'color', 'k', 'linew', 2);
        line([2 2], [mean_val_post+sd_val_post mean_val_post-sd_val_post], 'color','k','linew',1);
        line([1.95 2.05], [mean_val_post+sd_val_post mean_val_post+sd_val_post], 'color','k','linew',1);
        line([1.95 2.05], [mean_val_post-sd_val_post mean_val_post-sd_val_post], 'color','k','linew',1);
    end
    xlabel('Cohabitation Status');
    set(gca, 'xtick', [1 2], 'xticklabel', {'Before', 'After'});
    ylabel(fields{features_of_interest(i)});
    behav_name = fields{features_of_interest(i)};
    if contains(behav_name, 'Duration')
        ylim([0.04 0.12]);
        set(gca, 'ytick', linspace(0.04, 0.12, 5));
    elseif contains(behav_name, 'MedianFreq')
        ylim([30000 42000]); 
        set(gca, 'ytick', linspace(30000, 42000, 5));
    elseif contains(behav_name, 'HighFreq')
        ylim([35000 55000]); 
        set(gca, 'ytick', linspace(35000, 55000, 5));
    elseif contains(behav_name, 'LowFreq')
        ylim([25000 33000]); 
        set(gca, 'ytick', linspace(25000, 33000, 5));
    elseif contains(behav_name, 'Bandwidth')
        ylim([10000 22000]); 
        set(gca, 'ytick', linspace(10000, 22000, 5));
    elseif contains(behav_name, 'Slope')
        ylim([-20 40]); 
        set(gca, 'ytick', linspace(-20, 40, 5));
    elseif contains(behav_name, 'ModeNumHarmonics')
        ylim([0.1 0.5]); 
        set(gca, 'ytick', linspace(0.1, 0.5, 5));
    elseif contains(behav_name, 'NumSegments')
        ylim([1.4 2.2]); 
        set(gca, 'ytick', linspace(1.4, 2.2, 5));
    end
    set(gca, 'xlim', [.8 2.2]);
    ax = gca;
    ax.YColor = 'k';
    ax.XColor = 'k';
    ax.LineWidth = 1.0;
    saveas(gcf, fullfile(save_dir, sprintf('PreVsPostCohab-%s.jpg', behav_name)));
	saveas(gcf, fullfile(save_dir, sprintf('PreVsPostCohab-%s.svg', behav_name)));
end

%%Make a plot for number of USVs
figure('color','w'); hold on;
mean_vals = nanmean(num_usvs,1);
std_vals = nanstd(num_usvs,1)./sqrt(length(num_usvs));
if use_patch == 1
    patch([.7 .7 1.3 1.3], [mean_vals(1)-std_vals(1) mean_vals(1)+std_vals(1) mean_vals(1)+std_vals(1) mean_vals(1)-std_vals(1)], ...
        [253,245,230]./255);
    patch([1.7 1.7 2.3 2.3], [mean_vals(2)-std_vals(2) mean_vals(2)+std_vals(2) mean_vals(2)+std_vals(2) mean_vals(2)-std_vals(2)], ...
        [253,245,230]./255);
    line([.7 1.3], [mean_vals(1) mean_vals(1)], 'color','k','linew',3);
    line([1.7 2.3], [mean_vals(2) mean_vals(2)], 'color','k','linew',3);
else
    line([0.9 1.1], [mean_vals(1) mean_vals(1)], 'color', 'k', 'linew', 4);
    line([1 1], [mean_vals(1)-std_vals(1) mean_vals(1)+std_vals(1)], 'color','k','linew',2);
    line([.95 1.05], [mean_vals(1)-std_vals(1) mean_vals(1)-std_vals(1)], 'color','k','linew',2);
    line([.95 1.05], [mean_vals(1)+std_vals(1) mean_vals(1)+std_vals(1)], 'color','k','linew',2);
    line([1.9 2.1], [mean_vals(2) mean_vals(2)], 'color', 'k', 'linew', 4);
    line([2 2], [mean_vals(2)-std_vals(2) mean_vals(2)+std_vals(2)], 'color','k','linew',2);
    line([1.95 2.05], [mean_vals(2)-std_vals(2) mean_vals(2)-std_vals(2)], 'color','k','linew',2);
    line([1.95 2.05], [mean_vals(2)+std_vals(2) mean_vals(2)+std_vals(2)], 'color','k','linew',2);    
end
for j = 1:size(num_usvs,1)
    line([1 2], [num_usvs(j,1) num_usvs(j,2)], 'color','k', 'linew', 2);
    scatter([1 2]', [num_usvs(j,1) num_usvs(j,2)]', 60, repmat(scatter_colors(j,:),2,1), 'filled');
end
ylabel('Num USVs');
set(gca, 'xtick', [1 2], 'xticklabel', {'Before', 'After'});
xlabel('Cohabitation Status');
yvals = ylim;
ymin = yvals(1);
ymin = floor(ymin/1000)*1000;
ymax = ceil(yvals(2)/1000)*1000;
ylim([ymin ymax]);
set(gca, 'ytick', [ymin:1000:ymax]);

[h,p,ci,stats] = ttest(num_usvs(:,1), num_usvs(:,2))

saveas(gcf, fullfile(save_dir, 'NumVocsPreVsPostCohab.jpg'));
saveas(gcf, fullfile(save_dir, 'NumVocsPreVsPostCohab.svg'));
     
     

