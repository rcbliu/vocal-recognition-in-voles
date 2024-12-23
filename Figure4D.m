%%%script to plot violin plots of classifier accuracy when removing
%%%individual sinFM features from classifiers

clear
clc
close all

base_path = pwd;
save_dir = 'OutputPlots';

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

sin_fm_feats = {'Duration', 'AFM', 'FFM', 'OnsetFreq', 'Slope', 'Phi'};
ylabels = {'None', 'Dur', 'AFM', 'FFM', 'Onset', 'Slope', 'Phi'};

%%loads in AllConfus structure, which is recording x recording x shuffle x
%%feature (7x7x1000x6)
load('AllConfusMatrix_RemoveSinFMFeats');
s = load('ShuffConfusMatrix_ShuffSinFMFeats');
shuff_matrix = s.AllConfus;

%%Also get the non-shuffled data
s = load('ConfusMatrix_Day1');
class_matrix = s.AllConfus;
OverallAcc_class = zeros(1000,1);
for i = 1:1000
    acc_vals = 0;
    data_using = class_matrix(:,:,i);
    for k = 1:size(data_using,1)
        acc_vals = acc_vals+data_using(k,k);
    end
    acc_vals = acc_vals/size(data_using,1);
    OverallAcc_class(i,1) = acc_vals;
    clear acc_vals data_using
end

%%let's get the accuracy values

%%for each feature
OverallAcc = zeros(size(AllConfus,3),size(AllConfus,4));
for i = 1:size(AllConfus,4)
    data_tmp = AllConfus(:,:,:,i);
    %%for each shuffle
    for j = 1:size(AllConfus,3)
        data_using = data_tmp(:,:,j);
        %%for each recording
        acc_vals = 0;
        for k = 1:size(AllConfus,1)
            acc_vals = acc_vals+data_using(k,k);
        end
        acc_vals = acc_vals/size(data_using,1);
        OverallAcc(j,i) = acc_vals;
        clear acc_vals data_using
    end
    clear data_tmp
end

%%Add in the original classifier accuracy
OverallAcc = cat(2,OverallAcc_class,OverallAcc);

%%let's first plot the chance levels
figure('color','w'); hold on
%%make a patch for the shuffled results
shuff_acc = [];
for i = 1:size(shuff_matrix,3)
    shuff_data_tmp = shuff_matrix(:,:,i);
    shuff_acc_tmp = 0;
    for j = 1:size(shuff_matrix,2)
        shuff_acc_tmp = shuff_acc_tmp+shuff_data_tmp(j,j);
    end
    shuff_acc(i,1) = shuff_acc_tmp/size(shuff_matrix,2);
end
avg_shuff_acc = nanmean(shuff_acc);
std_shuff_acc = nanstd(shuff_acc);
P = patch('XData', [0 7 7 0], 'YData', ...
    [avg_shuff_acc-std_shuff_acc avg_shuff_acc-std_shuff_acc ...
    avg_shuff_acc+std_shuff_acc avg_shuff_acc+std_shuff_acc])
get(P)
P.FaceAlpha = 0.4;
P.FaceColor = 'r';
P.EdgeColor = 'none';
line([0 7], avg_shuff_acc+std_shuff_acc, 'color', 'r', 'linew', 1);
line([0 7], avg_shuff_acc-std_shuff_acc, 'color', 'r', 'linew', 1);

%%then plot our data
h1=distributionPlot_mrw_edits(gca,OverallAcc,'histOpt',1,'widthDiv',[2 1],'showMM',4,'addSpread',false,'distWidth',.7,'hist',1);
hold on
xvals = xlim;
line([xvals(1) xvals(2)], [1/7 1/7], 'color', 'r', 'linestyle', '--', 'linew', 2);

%%get stats for everything
for i = 1:size(OverallAcc,2)
    %%let's calculate a zscore
    avg_val = nanmean(OverallAcc(:,i));
    std_val = nanstd(OverallAcc(:,i));

    zscore(i,1) = (avg_val - nanmean(shuff_acc))/std_val;
end

p_vals = (1 - normcdf(zscore));
p_vals = cat(2,[1:length(p_vals)]', p_vals);
[B,I] = sort(p_vals(:,2), 'ascend');

p_val_labels = cat(2,p_vals(I,1),B);
p_val_labels(1,3) = 0.05;
for i = 2:length(p_val_labels)
    p_val_labels(i,3) = p_val_labels(i-1,3)/2;
end

set(gca, 'xticklabel', ylabels);
xtickangle(-45);

saveas(gcf, fullfile(save_dir, 'ClassAccuracy-RemoveSinFMFeats.jpg'));
saveas(gcf, fullfile(save_dir, 'ClassAccuracy-RemoveSinFMFeats.svg'));