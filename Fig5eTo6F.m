function [data_nosepoke, data_wc, nums_partner, nums_stranger, ...
    data_partner, data_stranger] = fn_quantify_side_times()

%%run stats, data for female mice attending to male-emitted sound
%%playback

clear
clc
close all

colors_pairs =[255 0 0; 0 204 0; 0 102 204; 204 0 204; 255 51 153; 255 128 0; 160 160 160]./255;
%%index the excel sheets
base_path = pwd;
floc = 'PlaybackBehavior';
save_dir = 'OutputPlots';
files = dir(fullfile(base_path, floc, '*.csv'));
fnames_tmp = {files.name}';

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

%%also pull in the data for partner location
rec_info = readtable('RecordingList.xlsx');

%%partner side for each recording
playback_sides = {'left',...
    'left',...
    'right',...
    'right',...
    'right',...
    'left'};


labels = {'Total number of occurences'	'Total duration (s)'	'Duration mean (s)'	'Duration std dev'}
Behaviors = {'Wall Touch left', 'Nose Poke left', 'Attend to right side', ...
    'Attend to left side', 'Wall Touch right', 'Nose Poke right'};

%%playback structures. 1 = USV, 0 = background.
A_structure = [1 0 0 1 0 1 1 0 0 1 0 1 0 1 1 0 1 0 1 0];
B_structure = [0 1 0 1 1 0 1 0 1 0 1 0 0 1 0 1 1 0 0 1];

head_to_head_times = find(A_structure == 1 & B_structure == 1);

times = [30:30:600];
times_A = times(A_structure==1)';
times_B = times(B_structure == 1)';
times_A = cat(2,times_A-30,times_A);
times_B = cat(2,times_B-30,times_B);
all_times(:,1) = [0:30:600-30];
all_times(:,2) = [30:30:600];

%%If this reads 'AB_LR', that means the left speaker is playing order A,
%%the right speaker is playing order B. 
playback_info = {'AB_LR'
    'BA_RL'
    'AB_LR'
    'AB_RL'
    'BA_RL'
    'AB_RL'};

%%let's also get times with only one speaker playing USVs, organized by
%%speaker
times_only_A = find(A_structure == 1 & B_structure == 0);
times_only_B = find(A_structure == 0 & B_structure == 1);
left_only_times = [];
right_only_times = [];

%%load in a file just to pull the behavior names
fname_tmp = fnames_tmp{1};
[path,fname,ext] = fileparts(fname_tmp);
[~,~,data] = xlsread(fullfile(base_path, floc, sprintf('%s.xlsx', fname)));
%%and naturally everything is stored in a super inconvenient format.
%%Because why not. Deal with that.
%%step 1; put all behavioral openings (all of the behaviors are state
%%behaviors) with their corresponding closes
data = data(17:end,:);
behav_names = unique(data(:,6));
num_behavs = zeros(length(playback_info),length(behav_names));
for i = 1:length(playback_info)
    playback_info_tmp = playback_info{i};
    if strcmp(playback_info_tmp(4), 'L')
        playback_left = playback_info_tmp(1);
        playback_right = playback_info_tmp(2);
    else
        playback_left = playback_info_tmp(2);
        playback_right = playback_info_tmp(1);
    end
    eval(sprintf('playback_times_left = times_%s;', playback_left));
    eval(sprintf('playback_times_right = times_%s;', playback_right));

    order_tmp = playback_info{i};
    if strcmp(order_tmp(4), 'L')
        if strcmp(order_tmp(1), 'A')
            left_only_times = times_only_A;
            right_only_times = times_only_B;
            playback_times_left = times_A;
            playback_times_right = times_B;
        else
            left_only_times = times_only_B;
            right_only_times = times_only_A;
            playback_times_left = times_B;
            playback_times_right = times_A;
        end
    elseif strcmp(order_tmp(4), 'R')
        if strcmp(order_tmp(1), 'A')
            right_only_times = times_only_A;
            left_only_times = times_only_B;
            playback_times_left = times_B;
            playback_times_right = times_A;
        else
            right_only_times = times_only_B;
            left_only_times = times_only_A;
            playback_times_left = times_A;
            playback_times_right = times_B;
        end
    end
    
    fname_tmp = fnames_tmp{i};
    [path,fname,ext] = fileparts(fname_tmp);
    [~,~,data] = xlsread(fullfile(floc, sprintf('%s.xlsx', fname)));
    %%and naturally everything is stored in a super inconvenient format.
    %%Because why not. Deal with that.
    %%step 1; put all behavioral openings (all of the behaviors are state
    %%behaviors) with their corresponding closes
    data = data(17:end,:);
    
    %%let's make a structure with the time of all nose pokes
    indices_nose_poke = strfind(data(:,6), 'Nose');
    indices_nose_poke = find(~cellfun(@isempty, indices_nose_poke));
    nose_poke_times_tmp = cell2mat(data(indices_nose_poke(1:2:end),1));
    nose_poke_times_tmp(:,2) = cell2mat(data(indices_nose_poke(2:2:end),1));
    AllPokeTimes{i} = nose_poke_times_tmp;
    clear indices_nose_poke nose_poke_times_tmp
    
    %%let's make a structure with the time of all wall touches
    indices_wall_touch = strfind(data(:,6), 'Wall');
    indices_wall_touch = find(~cellfun(@isempty, indices_wall_touch));
    wall_touch_times_tmp = cell2mat(data(indices_wall_touch(1:2:end),1));
    wall_touch_times_tmp(:,2) = cell2mat(data(indices_wall_touch(2:2:end),1));
    AllTouchTimes{i} = wall_touch_times_tmp;
    clear indices_wall_touch wall_touch_times_tmp
    
    for j = 3:length(behav_names)
        clear times_tmp
        data_tmp = data(strcmp(data(:,6), behav_names{j}),:);
        %%I am operating under the assumption that all consecutive pairs
        %%are start then stop
        times_tmp = data_tmp(1:2:end,1);
        times_tmp(:,2) = data_tmp(2:2:end,1);
        tot_time = [];
        if isempty(times_tmp)
            continue
        end
        times_tmp = cell2mat(times_tmp);
        
        %%built in self-check
        %%count the total duration of each behavior. Does it match the
        %%output from Boris?
        total_time(i,j) = nansum(times_tmp(:,2)-times_tmp(:,1));
        
        %%get rid of any times after 10 minutes
        idx_del_row = find(times_tmp(:,1) > 600);
        if ~isempty(idx_del_row)
            times_tmp(idx_del_row) = [];
        end
        idx_modify_time = find(times_tmp(:,2) > 600);
        if ~isempty(idx_modify_time)
            times_tmp(idx_modify_time,2) = 600;
        end
        %%calculate total time spent in this behavior when the relevant
        %%sound was on
        all_dur_tmp = 0;
        all_head_to_head_dur_left_tmp = 0;
        all_left_only_dur_tmp = 0;
        head_to_head_dur_tmp_left = 0;
        left_only_dur_tmp = 0;
        all_head_to_head_dur_right_tmp = 0;
        all_right_only_dur_tmp = 0;
        head_to_head_dur_tmp_right = 0;
        right_only_dur_tmp = 0;

        if contains(behav_names{j}, 'left')
            for k = 1:size(times_tmp,1)
                
                time_idx_all = find(times_tmp(k,1) > all_times(:,1) & ...
                    times_tmp(k,1) < all_times(:,2));
                time_idx = find(times_tmp(k,1) > playback_times_left(:,1) & ...
                    times_tmp(k,1) < playback_times_left(:,2));
                dur_tmp = 0;
                head_to_head_dur_tmp_left = 0;
                left_only_dur_tmp = 0;
                if ~isempty(time_idx)
                    %%If this behavior is split into two different time
                    %%blocks, only count the relevant block.
                    %%In theory, this becomes problematic if the next block
                    %%is also a sound block. Deal with that as well.
                    if times_tmp(k,2) > playback_times_left(time_idx,2) & ...
                            times_tmp(k,1) < playback_times_left(time_idx,2) & ...
                            playback_times_left(time_idx+1,1) ~= playback_times_left(time_idx,2)
                        dur_tmp = playback_times_left(time_idx,2) - times_tmp(k,1);
                    else
                        dur_tmp = times_tmp(k,2) - times_tmp(k,1);
                    end
                    
                    if ismember(time_idx_all, head_to_head_times)
                        if times_tmp(k,2) > playback_times_left(time_idx,2) & ...
                                times_tmp(k,1) < playback_times_left(time_idx,2) & ...
                                playback_times_left(time_idx+1,1) ~= playback_times_left(time_idx,2)
                            head_to_head_dur_tmp_left = playback_times_left(time_idx,2) - times_tmp(k,1);
                        else
                            head_to_head_dur_tmp_left = dur_tmp;
                        end
                    elseif ismember(time_idx_all, left_only_times)
                        if times_tmp(k,2) > playback_times_left(time_idx,2) & ...
                                times_tmp(k,1) < playback_times_left(time_idx,2) & ...
                                playback_times_left(time_idx+1,1) ~= playback_times_left(time_idx,2)
                            head_to_head_dur_tmp_left = playback_times_left(time_idx,2) - times_tmp(k,1);
                        else
                            left_only_dur_tmp = dur_tmp;
                        end
                    end
                end
                all_dur_tmp = all_dur_tmp + dur_tmp;
                all_head_to_head_dur_left_tmp = all_head_to_head_dur_left_tmp + ...
                    head_to_head_dur_tmp_left;
                all_left_only_dur_tmp = all_left_only_dur_tmp + left_only_dur_tmp;
                num_behavs(i,j) = num_behavs(i,j)+1;
                clear dur_tmp time_idx head_to_head_dur_tmp_left left_only_dur_tmp
            end
        else
            for k = 1:size(times_tmp,1)
                
                time_idx_all = find(times_tmp(k,1) > all_times(:,1) & ...
                    times_tmp(k,1) < all_times(:,2));
                time_idx = find(times_tmp(k,1) > playback_times_right(:,1) & ...
                    times_tmp(k,1) < playback_times_right(:,2));
                dur_tmp = 0;
                head_to_head_dur_tmp_right = 0;
                right_only_dur_tmp = 0;
                if ~isempty(time_idx)
                    if times_tmp(k,2) > playback_times_right(time_idx,2) & ...
                            times_tmp(k,1) < playback_times_right(time_idx,2) & ...
                            playback_times_right(time_idx+1,1) ~= playback_times_right(time_idx,2)
                        dur_tmp = playback_times_right(time_idx,2) - times_tmp(k,1);
                    else
                        dur_tmp = times_tmp(k,2) - times_tmp(k,1);
                    end
                    if ismember(time_idx_all, head_to_head_times)
                        if times_tmp(k,2) > playback_times_right(time_idx,2) & ...
                                times_tmp(k,1) < playback_times_right(time_idx,2) & ...
                                playback_times_right(time_idx+1,1) ~= playback_times_right(time_idx,2)
                            head_to_head_dur_tmp_right = playback_times_right(time_idx,2) - times_tmp(k,1);
                        else
                            head_to_head_dur_tmp_right = dur_tmp;
                        end
                    elseif ismember(time_idx_all, right_only_times)
                        if times_tmp(k,2) > playback_times_right(time_idx,2) & ...
                                times_tmp(k,1) < playback_times_right(time_idx,2) & ...
                                playback_times_right(time_idx+1,1) ~= playback_times_right(time_idx,2)
                            head_to_head_dur_tmp_right = playback_times_right(time_idx,2) - times_tmp(k,1);
                        else
                            right_only_dur_tmp = dur_tmp;
                        end
                    end
                end
                all_dur_tmp = all_dur_tmp + dur_tmp;
                all_head_to_head_dur_right_tmp = all_head_to_head_dur_right_tmp + ...
                    head_to_head_dur_tmp_right;
                all_right_only_dur_tmp = all_right_only_dur_tmp + right_only_dur_tmp;
                num_behavs(i,j) = num_behavs(i,j)+1;
                clear dur_tmp time_idx right_only_dur_tmp head_to_head_dur_tmp_right
            end
        end
        
        Head_to_head_dur_left_all(i,j) = all_head_to_head_dur_left_tmp;
        Head_to_head_dur_right_all(i,j) = all_head_to_head_dur_right_tmp;
        Right_only_dur_all(i,j) = all_right_only_dur_tmp;
        Left_only_dur_all(i,j) = all_left_only_dur_tmp;
        Durations(i,j) = all_dur_tmp;
        clear all_dur_tmp dur_tmp head_to_head_dur_tmp*
    end
end

%%let's plot the nose poke behavior over time
%%let's do a plot for individual lines and another with shaded error bars
figure('color','w');
%%let's break it into 1 minute bins
num_pokes = zeros(length(AllPokeTimes), 20);
dur_pokes = num_pokes;
num_wt = num_pokes;
dur_wt = dur_pokes;
for i = 1:length(AllPokeTimes)
    data_tmp = AllPokeTimes{i};
    data_wt = AllTouchTimes{i};
    count = 1;
    for j = 1:30:10*60
        idx_pokes = find(data_tmp(:,1) >= j & data_tmp(:,1) < j+30);
        if ~isempty(idx_pokes)
            num_pokes(i,count) = length(idx_pokes);
            dur_pokes(i,count) = nansum(data_tmp(idx_pokes,2) - data_tmp(idx_pokes,1));
        end
        idx_wt = find(data_wt(:,1) >= j & data_wt(:,1) < j+30);
        if ~isempty(idx_wt)
            num_wt(i,count) = length(idx_wt);
            dur_wt(i,count) = nansum(data_wt(idx_wt,2) - data_wt(idx_wt,1));
        end
        count = count+1;
    end
end

%%let's make bar plots
figure('color','w'); hold on
avg_dur_pokes = nanmean(dur_pokes,1)
err = nanstd(dur_pokes,1);
h = bar(avg_dur_pokes, 'facecolor', 'k', 'barw', 0.9);
h1 = errorbar([1:1:20], avg_dur_pokes, zeros(size(err,1),size(err,2)), err, '.'); % 'facecolor', 'k', 'barw', 1);
xlabel('Time (minutes)');
ylabel('Poke Duration (seconds)');
set(gca, 'xtick', [0:2:20], 'xticklabel', [0:1:10]);
%%add in errorbars
title_tmp = 'Bar Plot Duration Pokes Over Time';
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg', title_tmp)));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg', title_tmp)));


%%let's also do a bar plot for wall touch
figure('color','w'); hold on
avg_dur_touches = nanmean(dur_wt,1)
err = nanstd(dur_wt,1);
h = bar(avg_dur_touches, 'facecolor', 'k', 'barw', 0.9);
h1 = errorbar([1:1:20], avg_dur_touches, zeros(size(err,1),size(err,2)), err, '.'); % 'facecolor', 'k', 'barw', 1);
xlabel('Time (minutes)');
ylabel('Touch Duration (seconds)');
set(gca, 'xtick', [0:2:20], 'xticklabel', [0:1:10]);
%%add in errorbars
title_tmp = 'Bar Plot Duration WT Over Time';
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg', title_tmp)));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg', title_tmp)));

%%stats for wall touch
t2 = table([1:size(dur_wt,1)]',dur_wt(:,1),dur_wt(:,2),dur_wt(:,3),...
    dur_wt(:,4),dur_wt(:,5),dur_wt(:,6),dur_wt(:,7),...
    dur_wt(:,8),dur_wt(:,9),dur_wt(:,10),...
    dur_wt(:,11),dur_wt(:,12),dur_wt(:,13),...
    dur_wt(:,14),dur_wt(:,15),dur_wt(:,16),...
    dur_wt(:,17),dur_wt(:,18),dur_wt(:,19),dur_wt(:,20),...
    'VariableNames',{'id','bin1','bin2','bin3','bin4','bin5','bin6',...
    'bin7','bin8','bin9','bin10','bin11','bin12','bin13','bin14','bin15',...
    'bin16','bin17','bin18','bin19','bin20',});
times = table([0.5:0.5:10]','VariableNames',{'Bins'});

rm = fitrm(t2, 'bin1-bin20~id', 'WithinDesign', times);
ranovatbl = ranova(rm)
rm.Coefficients

%%let's run stats on the average duration of nose pokes over time
%%need a repeated measures ANOVA, so use ranovatbl = ranova(rm)
%%rm = fitrm(data,'output~predictor','WithinDesign',labels)
t1 = table([1:size(dur_pokes,1)]',dur_pokes(:,1),dur_pokes(:,2),dur_pokes(:,3),...
    dur_pokes(:,4),dur_pokes(:,5),dur_pokes(:,6),dur_pokes(:,7),...
    dur_pokes(:,8),dur_pokes(:,9),dur_pokes(:,10),...
    dur_pokes(:,11),dur_pokes(:,12),dur_pokes(:,13),...
    dur_pokes(:,14),dur_pokes(:,15),dur_pokes(:,16),...
    dur_pokes(:,17),dur_pokes(:,18),dur_pokes(:,19),dur_pokes(:,20),...
    'VariableNames',{'id','bin1','bin2','bin3','bin4','bin5','bin6',...
    'bin7','bin8','bin9','bin10','bin11','bin12','bin13','bin14','bin15',...
    'bin16','bin17','bin18','bin19','bin20',});
times = table([0.5:0.5:10]','VariableNames',{'Bins'});

rm = fitrm(t1, 'bin1-bin20~id', 'WithinDesign', times);
ranovatbl = ranova(rm)

%%now sort the data back out based on mate vs other side
indices_left = find(contains(behav_names,'left'));
indices_right = find(contains(behav_names,'right'));

%%behavior structure = number of attends, num nose pokes, num WC
for i = 1:length(playback_sides)
    if strcmp(playback_sides{i}, 'left')
        data_partner(i,:) = Durations(i,indices_left);
        data_stranger(i,:) = Durations(i,indices_right);
        nums_partner(i,:) = num_behavs(i,indices_left);
        nums_stranger(i,:) = num_behavs(i,indices_right);
        head_to_head_data_partner(i,:) = Head_to_head_dur_left_all(i,indices_left);
        head_to_head_data_stranger(i,:) = Head_to_head_dur_right_all(i,indices_right);
        partner_only_playback(i,:) = Left_only_dur_all(i,indices_left);
        stranger_only_playback(i,:) = Right_only_dur_all(i,indices_right);
    else
        data_partner(i,:) = Durations(i,indices_right);
        data_stranger(i,:) = Durations(i,indices_left);
        nums_partner(i,:) = num_behavs(i,indices_right);
        nums_stranger(i,:) = num_behavs(i,indices_left);
        head_to_head_data_partner(i,:) = Head_to_head_dur_right_all(i,indices_right);
        head_to_head_data_stranger(i,:) = Head_to_head_dur_left_all(i,indices_left);
        partner_only_playback(i,:) = Right_only_dur_all(i,indices_right);
        stranger_only_playback(i,:) = Left_only_dur_all(i,indices_left);
    end
end
%%head to head files have each row being a recording, each column being a
%%behavior (attend, nose poke, wall touch)
nums_partner(end,2) = 1;
data_partner(end,2) = 0.5;

figure('color','w'); hold on
colors = distinguishable_colors(length(playback_sides),'w');
h1 = scatter([1:1:6], data_partner(:,1), 30, colors, 'filled');
set(h1, 'markeredgecolor', 'k', 'linew', 0.5);
%%make lighter colors
cmap = colors;
newcolors = brighten(cmap,0.7);
%%plot stranger data
h2 = scatter([1:1:6], data_stranger(:,1), 30, newcolors, 'filled');

%%Plot proportion of nose pokes on partner side. 
%%Nose pokes are in column 2, followed by wall touches
data_plot = data_partner(:,2) ./ (data_partner(:,2)+data_stranger(:,2));
figure('color','w'); bar(data_plot);
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2);

%%stats nose poke time
[h,p,ci,stats] = ttest2(data_plot, 0.5)

data_plot = data_partner(:,2) ./ (data_partner(:,2)+data_stranger(:,2));
idx_zero = find(data_plot == 0 | isnan(data_plot));
%%if neither animal poked, do we refer to it as 0 or .5 or 1? I am voting
%%0.5
for i = 1:length(idx_zero)
    if data_partner(idx_zero(i),2) + data_stranger(idx_zero(i),2) == 0
        data_plot(idx_zero(i)) = 0.5;
    end
end
data_plot(:,2) = 1-data_plot;
data_plot(end,:) = [0.5 0.5];
figure('color','w'); h = bar(data_plot, 'stacked');
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2);
title('Proportion of nose poke time on partner side during sound playback');
box off

hold on
range_partner = prctile(data_plot(:,1), [25 50 75]);
range_stranger = prctile(data_plot(:,2), [25 50 75]);

colors_bars = get(h,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));
% d = computeCohen_d(data_partner(:,2)./(data_partner(:,2)+data_stranger(:,2)), ...
%     repmat(0.5, length(data_partner(:,2)),1))

%%try stats wall touch time
[h,p,ci,stats] = ttest(data_partner(:,3)./(data_partner(:,3)+data_stranger(:,3)), 0.5)

data_plot = data_partner(:,3) ./ (data_partner(:,3)+data_stranger(:,3));
figure('color','w'); bar(data_plot);
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2);
title('Proportion of wall touch time on partner side during sound playback');

%%try stats nose poke time
[h,p,ci,stats] = ttest(data_plot,0.5)

%%look at proportion of number of times animal nose pokes on partner side
data_plot = nums_partner(:,2)./(nums_partner(:,2)+nums_stranger(:,2));
idx_zero = find(data_plot == 0 | isnan(data_plot));
data_plot(idx_zero) = 0.5;
figure('color','w'); bar(data_plot);
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2);
title('Proportion of num nose pokes to partner side during sound playback');
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

%%try stats
[h,p,ci,stats] = ttest(nums_partner(:,2)./(nums_partner(:,2)+nums_stranger(:,2)), 0.5)

%%also do total num
figure('color','w'); h = bar(cat(2,nums_partner(:,2), nums_stranger(:,2)), 'stacked');
xvals = xlim;
title_tmp = 'Num nose pokes to partner side during sound playback';
title(title_tmp);
box off
legend({'Partner', 'Stranger'});
ylabel('Number Nose Pokes')
xlabel('Pair ID')

hold on
range_partner = prctile(nums_partner(:,2), [25 50 75]);
range_stranger = prctile(nums_stranger(:,2), [25 50 75]);

colors_bars = get(h,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));


[h,p,stats] = ttest(nums_partner(:,2), nums_stranger(:,2));

%%also do total dur
figure('color','w'); h = bar(cat(2,data_partner(:,2), data_stranger(:,2)), 'stacked');
xvals = xlim;
title_tmp = 'Dur nose pokes to partner side during sound playback';
title(title_tmp);
box off
ylabel('Duration Nose Pokes')
xlabel('Pair ID')
legend({'Partner', 'Stranger'});

hold on
range_partner = prctile(data_partner(:,2), [25 50 75]);
range_stranger = prctile(data_stranger(:,2), [25 50 75]);

colors_bars = get(h,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

[h,p,ci,stats] = ttest(data_partner(:,2), data_stranger(:,2));

%%look at proportion of number of times animal wall climbs on partner side
data_plot = nums_partner(:,3)./(nums_partner(:,3)+nums_stranger(:,3));
figure('color','w'); bar(data_plot);
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2);
title('Proportion of wall touches to partner side during sound playback');
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

%%try stats
[h,p] = ttest(0.5 - nums_partner(:,3)./(nums_partner(:,3)+nums_stranger(:,3)))


%%make new plot. Dot plot for the poster with both noste pokes and wall
%%touches
data_wc = data_partner(:,3)./(data_partner(:,3)+data_stranger(:,3));
data_nosepoke = data_partner(:,2)./(data_partner(:,2)+data_stranger(:,2));
idx_zero = find(data_nosepoke == 0 | isnan(data_nosepoke));

for i = 1:length(idx_zero)
    if data_partner(idx_zero(i),2) + data_stranger(idx_zero(i),2) == 0
        data_nosepoke(idx_zero(i)) = 0.5;
    end
end
figure('color','w'); hold on
%%let's plot mean and std as boxes
mean_wc = nanmean(data_wc); 
std_wc = nanstd(data_wc)./sqrt(length(data_wc));
mean_np = nanmean(data_nosepoke);
std_np = nanstd(data_nosepoke)./sqrt(length(data_nosepoke));
% patch([0.8 0.8 1.2 1.2], [mean_wc+std_wc mean_wc-std_wc mean_wc-std_wc mean_wc+std_wc], ...
%     [255,245,238]./255);
% patch([1.8 1.8 2.2 2.2], [mean_np+std_np mean_np-std_np mean_np-std_np mean_np+std_np], ...
%     [255,245,238]./255);
line([1 1], [mean_wc+std_wc mean_wc-std_wc], 'color','k','linew',2);
line([.95 1.05], [mean_wc+std_wc mean_wc+std_wc], 'color','k','linew',2);
line([.95 1.05], [mean_wc-std_wc mean_wc-std_wc], 'color','k','linew',2);

line([1.5 1.5], [mean_np+std_np mean_np-std_np], 'color','k','linew',2);
line([1.45 1.55], [mean_np+std_np mean_np+std_np], 'color','k','linew',2);
line([1.45 1.55], [mean_np-std_np mean_np-std_np], 'color','k','linew',2);

line([.9 1.1], [mean_wc mean_wc], 'color','k','linew',4);
line([1.4 1.6], [mean_np mean_np], 'color','k','linew',4);
    
data_nosepoke(end) = 0.5;
line([0.7 2.3], [0.5 0.5], 'color', 'r', 'linew', 3);
for i = 1:length(nums_partner)
    line([1 1.5], [data_wc(i) data_nosepoke(i)], 'color', [200 200 200]./255, 'linew', 2);
    scatter([1 1.5], [data_wc(i) data_nosepoke(i)], 60, repmat(colors_pairs(i,:),2,1), 'filled');
end
set(gca, 'xlim', [0.7 1.8], 'ylim', [0 1]);
set(gca, 'ytick', [0:.25:1]);
set(gca, 'xtick', [1 1.5], 'xticklabel', {'WallTouch', 'NosePoke'});
ylabel('Proportion of Time Towards Partner')
title('Proportion of Behavior Time Directed Toward Partner - DotPlot');
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string')))); 


%%Another analysis: what was the poking false alarm rate? AKA, how many
%%times/how much time did the animal nose poke towards no sound.
%%behaviors are rows; labels are columns
%%This is coordinated data from excel behav counts
clear data_tmp
data_tmp{1,1} = [7	42.777	6.111	6.9
    18	52.953	2.942	2.429
    29	68.471	2.361	2.644
    15	96.232	6.415	8.761
    9	12.413	1.379	0.492
    8	3.271	0.409	0.282];

data_tmp{2,1} = [4	10.644	2.661	2.484
    8	11.344	1.418	1.261
    31	88.156	2.844	4.149
    24	44.074	1.836	2.266
    11	39.707	3.61	3.149
    8	11.01	1.376	1.862]

data_tmp{3,1} = [9	15.782	1.754	1.573
    9	18.954	2.106	3.231
    31	168.47	5.435	8.235
    22	62.596	2.845	3.519
    11	55.988	5.09	4.922
    17	84.643	4.979	4.78
    ];

data_tmp{4,1} = [1	4.137	4.137	0
    4	1.068	0.267	0.061
    27	49.385	1.829	2.457
    15	25.859	1.724	1.774
    6	16.684	2.781	2.946
    12	8.376	0.698	0.538
    ];

data_tmp{5,1} = [1	0.667	0.667	0
    0	0	0	0
    10	43.243	4.324	5.298
    7	5.772	0.825	0.379
    2	2.37	1.185	1.581
    6	13.714	2.286	2.352
    ];

data_tmp{6,1} = [1	0.767	0.767	0
    0	0	0	0
    30	45.48	1.516	1.373
    31	48.981	1.58	1.893
    3	4.972	1.657	1.38
    1	0.167	0.167	0
    ];
%%nose poke left = 2, right = 6
for i = 1:length(data_tmp)
    partner_side_tmp = playback_sides{i};
    if strcmp(partner_side_tmp, 'left')
        side_times_nose(i,1) = data_tmp{i,1}(2,2);
        side_times_nose(i,2) = data_tmp{i,1}(6,2);
    elseif strcmp(partner_side_tmp, 'right')
        side_times_nose(i,1) = data_tmp{i,1}(6,2);
        side_times_nose(i,2) = data_tmp{i,1}(2,2);
    end
end
for i = 1:length(data_tmp)
    partner_side_tmp = playback_sides{i};
    if strcmp(partner_side_tmp, 'left')
        side_times_touch(i,1) = data_tmp{i,1}(1,2);
        side_times_touch(i,2) = data_tmp{i,1}(5,2);
    elseif strcmp(partner_side_tmp, 'right')
        side_times_touch(i,1) = data_tmp{i,1}(5,2);
        side_times_touch(i,2) = data_tmp{i,1}(1,2);
    end
end


time_poking_no_sound(:,1) = side_times_nose(:,1)-data_partner(:,2);
time_poking_no_sound(:,2) = side_times_nose(:,2)-data_stranger(:,2);
time_poking_no_sound(time_poking_no_sound<0) = 0; %%get rid of rounding errors

time_wc_no_sound(:,1) = side_times_touch(:,1) - data_partner(:,3);
time_wc_no_sound(:,2) = side_times_touch(:,2) - data_stranger(:,3);
time_wc_no_sound(time_wc_no_sound<0) = 0;

%%make plot of proportion of time per side that is a false alarm (partner
%%left, stranger right)
clear false_alarm_prop
false_alarm_time = nansum(time_poking_no_sound,2);
false_alarm_prop = false_alarm_time./nansum(side_times_nose,2);
false_alarm_prop(end) = 0.5; 
false_alarm_prop = cat(2,1-false_alarm_prop,false_alarm_prop); %%put true pokes before false pokes

figure('color','w'); ba = bar(false_alarm_prop, 'stacked');
ba(1).FaceColor = [107 142 35]/255;
ba(2).FaceColor = [1 1 1]*0.8;
legend({'TruePositive', 'FalseNegative'});
% set(gca, 'xlim', [0.5 6.5]);
xlabel('Pair ID');
ylabel('Proportion');
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2)

hold on
range_partner = prctile(false_alarm_prop(:,1), [25 50 75]);
range_stranger = prctile(false_alarm_prop(:,2), [25 50 75]);

colors_bars = get(ba,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

[h,p,ci,stats] = ttest(false_alarm_prop(:,1), 0.5);

saveas(gcf, fullfile(save_dir, 'PropTrueVsFalsePositiveNP.jpg'));
saveas(gcf, fullfile(save_dir, 'PropTrueVsFalsePositiveNP.svg'));

%%let's also do wall touch
clear false_alarm_prop
false_alarm_time = nansum(time_wc_no_sound,2);
false_alarm_prop = false_alarm_time./nansum(side_times_touch,2);
false_alarm_prop(end) = 0.5; 
false_alarm_prop = cat(2,1-false_alarm_prop,false_alarm_prop); %%put true pokes before false pokes

figure('color','w'); ba = bar(false_alarm_prop, 'stacked');
ba(1).FaceColor = [107 142 35]/255;
ba(2).FaceColor = [1 1 1]*0.8;
legend({'TruePositive', 'FalseNegative'});
% set(gca, 'xlim', [0.5 6.5]);
xlabel('Pair ID');
ylabel('Proportion');
xvals = xlim;
line([xvals(1) xvals(2)], [0.5 0.5], 'color', 'k', 'linew', 2)

hold on
range_partner = prctile(false_alarm_prop(:,1), [25 50 75]);
range_stranger = prctile(false_alarm_prop(:,2), [25 50 75]);

colors_bars = get(ba,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

[h,p,ci,stats] = ttest(false_alarm_prop(:,1), 0.5);

saveas(gcf, fullfile(save_dir, 'PropTrueVsFalsePositiveWallTouch.jpg'));
saveas(gcf, fullfile(save_dir, 'PropTrueVsFalsePositiveWallTouch.svg'));

%%Also add in quantification of what happened at times when both speakers
%%were playing USVs: in a head-to-head competition, did the animal seem to
%%prefer her own male's USVs?

prop_head_to_head_time_partner = head_to_head_data_partner(:,2)./...
    (head_to_head_data_partner(:,2)+head_to_head_data_stranger(:,2));
prop_head_to_head_time_stranger = head_to_head_data_stranger(:,2)./...
    (head_to_head_data_partner(:,2)+head_to_head_data_stranger(:,2));

prop_head_to_head_time_partner(isnan(prop_head_to_head_time_partner)) = 0.5;
prop_head_to_head_time_stranger(isnan(prop_head_to_head_time_stranger)) = 0.5;

figure('color','w'); ba = bar(cat(2,prop_head_to_head_time_partner,prop_head_to_head_time_stranger), 'stacked');
legend({'Partner', 'Stranger'});
ylabel('Prop Head-to-head Time');
xlabel('RecID');
xvals = xlim;
line([xvals(1) xvals(2)], [.5 .5], 'color', 'k', 'linew', 2);

hold on
range_partner = prctile(prop_head_to_head_time_partner, [25 50 75]);
range_stranger = prctile(prop_head_to_head_time_stranger, [25 50 75]);

colors_bars = get(ba,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

saveas(gcf, fullfile(save_dir, 'PropHeadToHeadTimeTowardsPartner.jpg'));
saveas(gcf, fullfile(save_dir, 'PropHeadToHeadTimeTowardsPartner.svg'));

[h,p,ci,stats] = ttest(prop_head_to_head_time_partner, 0.5)

prop_only_one_speaker_partner = partner_only_playback(:,2)./...
    (partner_only_playback(:,2)+stranger_only_playback(:,2));
prop_only_one_speaker_stranger = stranger_only_playback(:,2)./...
    (partner_only_playback(:,2)+stranger_only_playback(:,2));
prop_only_one_speaker_partner(isnan(prop_only_one_speaker_partner)) = 0.5;
prop_only_one_speaker_stranger(isnan(prop_only_one_speaker_stranger)) = 0.5;

figure('color','w'); ba = bar(cat(2,prop_only_one_speaker_partner,prop_only_one_speaker_stranger), 'stacked');
legend({'Partner Only', 'Stranger Only'}, 'location', 'best');
box off
title('Proportion of NosePokeTime - Only1Speaker Playback')
ylabel('Proportion of Time');
xlabel('Pair ID')
xvals = xlim;
line([xvals(1) xvals(2)], [.5 .5], 'color', 'k', 'linew', 2);

hold on
range_partner = prctile(prop_only_one_speaker_partner, [25 50 75]);
range_stranger = prctile(prop_only_one_speaker_stranger, [25 50 75]);

colors_bars = get(ba,'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

[h,p,ci,stats] = ttest(prop_only_one_speaker_partner, 0.5)

%%let's also add a stacked barplot showing all USV-directed behavior during
%%head-to-head and during solo USVs
idx_add = find(head_to_head_data_partner(:,2) == 0);
head_to_head_data_partner(idx_add,2) = 0.5;
figure('color','w'); hold on
one_speaker_only = partner_only_playback(:,2)+stranger_only_playback(:,2);
two_speakers = head_to_head_data_partner(:,2)+head_to_head_data_stranger(:,2);
ba = bar(cat(2,one_speaker_only,two_speakers), 'stacked');
legend({'PartnerOnly', 'HeadToHead'}, 'location', 'best');
title('TotalDurNosePokes-HeadToHeadVsSolo')
set(gca, 'ylim', [0 70]);

hold on
range_partner = prctile(one_speaker_only, [25 50 75]);
range_stranger = prctile(two_speakers, [25 50 75]);

colors_bars = get(ba,'FaceColor');
ba(1).FaceColor = [73 111 255]./255;
ba(2).FaceColor = [153 166 255]./255;

colors_bars = get(ba, 'FaceColor');

h = patch('XData', [6.5 6.9 6.9 6.5], 'YData', [range_partner(1) range_partner(1) range_partner(3) range_partner(3)]);
h1 = patch('XData', [7.1 7.5 7.5 7.1], 'YData', [range_stranger(1) range_stranger(1) range_stranger(3) range_stranger(3)]);
h.FaceColor = colors_bars{1};
h1.FaceColor = colors_bars{2};

line([6.5 6.9], [range_partner(2) range_partner(2)], 'color', 'k', 'linew', 2);
line([7.1 7.5], [range_stranger(2) range_stranger(2)], 'color', 'k', 'linew', 2);

[h,p,ci,stats] = ttest(one_speaker_only,two_speakers)

text(3, 50, sprintf('p=%.2f', p));

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

%%let's also add a stacked barplot showing partner-directed behavior during
%%head-to-head and during solo USVs
partner_only_playback(end,2) = 1;
idx_add = find(head_to_head_data_partner(:,2) == 0);
head_to_head_data_partner(idx_add,2) = 0.5;
figure('color','w'); hold on
h = bar(cat(2,partner_only_playback(:,2),head_to_head_data_partner(:,2)), 'stacked');
legend({'PartnerOnly', 'HeadToHead'}, 'location', 'best');
title('TotalDurNosePokes-HeadToHeadVsSolo')
set(gca, 'ylim', [0 70]);


[h,p,ci,stats] = ttest(partner_only_playback(:,2),head_to_head_data_partner(:,2))

text(3, 50, sprintf('p=%.2f', p));

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));

clear partner_only_prop head_to_head_prop
figure('color','w'); hold on
partner_only_prop = partner_only_playback(:,2)./(partner_only_playback(:,2)+head_to_head_data_partner(:,2));
head_to_head_prop = head_to_head_data_partner(:,2)./(partner_only_playback(:,2)+head_to_head_data_partner(:,2));
h = bar(cat(2,partner_only_prop(:,1),head_to_head_prop(:,1)), 'stacked');
legend({'ParnterOnly', 'HeadToHead'}, 'location', 'best');
title('Prop Time Partner Only vs HeadToHead');
xvals = xlim;
line([xvals(1) xvals(2)], [.5 .5], 'color', 'k', 'linew', 2);

[h,p,ci,stats] = ttest(partner_only_prop,0.5);

saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, ...
    sprintf('%s.svg',get(get(gca,'title'),'string'))));
end

