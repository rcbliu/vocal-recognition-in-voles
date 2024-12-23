%%%Want to make a figure showing the cdf plots for pretest to prepair,
%%%pretest+prepair, and prepair to postpair classification. %
%%%As such, copying over a bunch of other code just to get all the results
%%%in one space for plotting purposes.

clear
clc
close all

%%script to use prepair data to predict pretest identity using
%%sinFM features
rng default

base_path = pwd;
sinfm_loc = 'tsne_Data';
save_dir = 'OutputPlots';

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

sin_fm_feats = {'Duration', 'AFM', 'FFM', 'OnsetFreq', 'Slope', 'Phi'};

pairings = {'618-620', '626'; ...
    '622', '627'; ...
    '624', '628'; ...
    '625', '629'; ...
    '656', '642'; ...
    '658', '643'; ...
    '659', '644'};


AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);

min_num_vocs = 125;
testing_size = floor(min_num_vocs*.25);
training_size = floor(min_num_vocs*.75);

tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    %%pretest is line 1, prepair is line 2
    for j = 1:size(pairings,1)
        clear data_poss
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,1})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,1})));
        end
        data_train = s.param_mx(:,1:6);
        for k = 1:size(data_train,2)
            data_train(:,k) = mat2gray(data_train(:,k));
        end
        clear s 

        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,2})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,2})));
        end
        data_test = s.param_mx(:,1:6);
        for k = 1:size(data_test,2)
            data_test(:,k) = mat2gray(data_test(:,k));
        end
        clear s 

        clear randvals
        randvals_train = randperm(length(data_train),min_num_vocs);
        randvals_test = randperm(length(data_test),min_num_vocs);
        training_data = cat(1,training_data,data_train(randvals_train(1:training_size),:));
        testing_data = cat(1,testing_data,data_test(randvals_test(1:testing_size),:));
        individual_ids_training(end+1:end+training_size,1) = repmat(j,training_size,1);
        individual_ids_testing(end+1:end+testing_size,1) = repmat(j,testing_size,1);
    end
    Mdl = fitcecoc(training_data,individual_ids_training);
    [label,score] = predict(Mdl,testing_data);
    prop_corr(i,1) = length(find(individual_ids_testing == label))./length(label);
    AllLabels{i,1} = cat(2,label,individual_ids_testing);
    confus_matrix = zeros(7,7);
    for j = 1:length(AllLabels{i,1})
        confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1)) = confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1))+1;
    end
    confus_matrix = bsxfun(@rdivide, confus_matrix, sum(confus_matrix,2));
    AllConfus(:,:,i) = confus_matrix;

    clear randvals randvals_* data_train data_test labels_train labels_test
    clear Mdl label score confus_matrix

    if rem(i,10) == 0
        fprintf('%2.0f%% complete %.2f\n', (i/1000)*100, toc)
    end
end
toc

class_acc = zeros(7,1000);
for i = 1:length(AllConfus)
    for j = 1:7
        class_acc(j,i) = AllConfus(j,j,i);
    end
end
class_acc = class_acc.*1000;
overall_class_acc = nansum(class_acc,1)/7000;

class_acc_pretest_to_prepair = class_acc;
overall_class_acc_pretest_to_prepair = overall_class_acc;


%%Step 1: let's do the same for the prepair and pretest data combined
rng default

AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);

min_num_vocs = 125;
testing_size = floor(min_num_vocs*.25);
training_size = floor(min_num_vocs*.75);

tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    %%pretest is line 1, prepair is line 2
    for j = 1:size(pairings,1)
        %%pretest data
        clear data_poss
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,1})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,1})));
        end
        data_train = s.param_mx(:,1:6);
        for k = 1:size(data_train,2)
            data_train(:,k) = mat2gray(data_train(:,k));
        end
        clear s 

        %%prepair data
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,2})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,2})));
        end
        data_train2 = s.param_mx(:,1:6);
        for k = 1:size(data_train2,2)
            data_train2(:,k) = mat2gray(data_train2(:,k));
        end
        clear s 



        clear randvals
        randvals_train = randperm(length(data_train),min_num_vocs);
        randvals_train2 = randperm(length(data_train2),min_num_vocs);

        training_data = cat(1,training_data,data_train(randvals_train(1:floor(training_size/2)),:),data_train2(randvals_train(1:floor(training_size/2)),:));
        testing_data = cat(1,testing_data,...
            data_train(randvals_train(ceil(training_size/2)+1:ceil(training_size/2)+floor(testing_size/2)),:),...
            data_train2(randvals_train(ceil(training_size/2)+1:ceil(training_size/2)+floor(testing_size/2)),:));
        individual_ids_training(end+1:end+floor(training_size/2)*2,1) = repmat(j,floor(training_size/2)*2,1);
        individual_ids_testing(end+1:end+floor(testing_size/2)*2,1) = repmat(j,floor(testing_size/2)*2,1);
    end
    Mdl = fitcecoc(training_data,individual_ids_training);
    [label,score] = predict(Mdl,testing_data);
    prop_corr(i,1) = length(find(individual_ids_testing == label))./length(label);
    AllLabels{i,1} = cat(2,label,individual_ids_testing);
    confus_matrix = zeros(7,7);
    for j = 1:length(AllLabels{i,1})
        confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1)) = confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1))+1;
    end
    confus_matrix = bsxfun(@rdivide, confus_matrix, sum(confus_matrix,2));
    AllConfus(:,:,i) = confus_matrix;

    clear randvals randvals_* data_train data_test labels_train labels_test
    clear Mdl label score confus_matrix

    if rem(i,10) == 0
        fprintf('%2.0f%% complete %.2f\n', (i/1000)*100, toc)
    end
end
toc

class_acc = zeros(7,1000);
for i = 1:length(AllConfus)
    for j = 1:7
        class_acc(j,i) = AllConfus(j,j,i);
    end
end
class_acc = class_acc.*1000;
overall_class_acc = nansum(class_acc,1)/7000;

class_acc_pretest_and_prepair = class_acc;
overall_class_acc_pretest_and_prepair = overall_class_acc;

%% 
%%3) Use pre data to predict post identity
rng default

AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);

min_num_vocs = 125;
testing_size = floor(min_num_vocs*.25);
training_size = floor(min_num_vocs*.75);

tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    %%pretest is line 1, prepair is line 2
    for j = 1:size(pairings,1)
        %%pretest data
        clear data_poss
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,1})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,1})));
        end
        data_train = s.param_mx(:,1:6);
        for k = 1:size(data_train,2)
            data_train(:,k) = mat2gray(data_train(:,k));
        end
        clear s 

        %%prepair data
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,2})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,2})));
        end
        data_train2 = s.param_mx(:,1:6);
        for k = 1:size(data_train2,2)
            data_train2(:,k) = mat2gray(data_train2(:,k));
        end
        clear s 

        %%postpair data
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,2})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,2})));
        end
        data_test = s.param_mx(:,1:6);
        for k = 1:size(data_test,2)
            data_test(:,k) = mat2gray(data_test(:,k));
        end
        clear s 

        clear randvals
        randvals_train = randperm(length(data_train),min_num_vocs);
        randvals_train2 = randperm(length(data_train2),min_num_vocs);
        randvals_test = randperm(length(data_test),min_num_vocs);
        training_data = cat(1,training_data,data_train(randvals_train(1:floor(training_size/2)),:),data_train2(randvals_train(1:floor(training_size/2)),:));
        testing_data = cat(1,testing_data,data_test(randvals_test(1:testing_size),:));
        individual_ids_training(end+1:end+floor(training_size/2)*2,1) = repmat(j,floor(training_size/2)*2,1);
        individual_ids_testing(end+1:end+testing_size,1) = repmat(j,testing_size,1);
    end
    Mdl = fitcecoc(training_data,individual_ids_training);
    [label,score] = predict(Mdl,testing_data);
    prop_corr(i,1) = length(find(individual_ids_testing == label))./length(label);
    AllLabels{i,1} = cat(2,label,individual_ids_testing);
    confus_matrix = zeros(7,7);
    for j = 1:length(AllLabels{i,1})
        confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1)) = confus_matrix(AllLabels{i,1}(j,2),AllLabels{i,1}(j,1))+1;
    end
    confus_matrix = bsxfun(@rdivide, confus_matrix, sum(confus_matrix,2));
    AllConfus(:,:,i) = confus_matrix;

    clear randvals randvals_* data_train data_test labels_train labels_test
    clear Mdl label score confus_matrix

    if rem(i,10) == 0
        fprintf('%2.0f%% complete %.2f\n', (i/1000)*100, toc)
    end
end
toc

class_acc = zeros(7,1000);
for i = 1:length(AllConfus)
    for j = 1:7
        class_acc(j,i) = AllConfus(j,j,i);
    end
end
class_acc = class_acc.*1000;
overall_class_acc = nansum(class_acc,1)/7000;

class_acc_pre_to_post = class_acc;
overall_class_acc_pre_to_post = overall_class_acc;

%%do a shuffle classifier as well for comparison
tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    %%pretest is line 1, prepair is line 2
    for j = 1:size(pairings,1)
        clear data_poss
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,1})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,1})));
        end
        data_train = s.param_mx(:,1:6);
        for k = 1:size(data_train,2)
            data_train(:,k) = mat2gray(data_train(:,k));
        end
        clear s 
        
        try
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_ch1T0000%s.mat', pairings{j,2})));
        catch
            s = load(fullfile(sinfm_loc, sprintf('all_call_data_T0000%s.mat', pairings{j,2})));
        end
        data_test = s.param_mx(:,1:6);
        for k = 1:size(data_test,2)
            data_test(:,k) = mat2gray(data_test(:,k));
        end
        clear s 

        
        clear randvals
        randvals_train = randperm(length(data_train),min_num_vocs);
        randvals_test = randperm(length(data_test), min_num_vocs);
        data_using = cat(1,data_train(randvals_train,:),data_test(randvals_test,:));

        randvals = randperm(length(data_using),length(data_using));

        data_using = data_using(randvals,:);
        training_data = cat(1,training_data,data_using(1:training_size,:));
        testing_data = cat(1,testing_data,data_using(training_size+1:training_size+testing_size,:));
        individual_ids_training(end+1:end+training_size,1) = repmat(j,training_size,1);
        individual_ids_testing(end+1:end+testing_size,1) = repmat(j,testing_size,1);
    end
    individual_ids_testing = individual_ids_testing(randperm(length(individual_ids_testing)));
    Mdl = fitcecoc(training_data,individual_ids_training);
    [label,score] = predict(Mdl,testing_data);
    prop_corr_shuff(i,1) = length(find(individual_ids_testing == label))./length(label);
    AllLabels_shuff{i,1} = cat(2,label,individual_ids_testing);
    confus_matrix = zeros(7,7);
    for j = 1:length(AllLabels{i,1})
        confus_matrix(AllLabels_shuff{i,1}(j,2),AllLabels_shuff{i,1}(j,1)) = confus_matrix(AllLabels_shuff{i,1}(j,2),AllLabels{i,1}(j,1))+1;
    end
    confus_matrix = bsxfun(@rdivide, confus_matrix, sum(confus_matrix,2));
    AllConfus_shuff(:,:,i) = confus_matrix;
    
    clear randvals randvals_* data_train data_test labels_train labels_test
    clear Mdl label score confus_matrix
        
    if rem(i,10) == 0
        fprintf('%2.0f%% complete %.2f\n', (i/1000)*100, toc)
    end
end
toc


%%let's get the accuracy.
shuff_acc = zeros(7,1000);
for i = 1:length(AllConfus_shuff)
    for j = 1:7
        shuff_acc(j,i) = AllConfus_shuff(j,j,i);
    end
end
shuff_acc = shuff_acc.*1000;
overall_shuff_acc = nansum(shuff_acc,1)/7000;

figure('color','w'); hold on
h = cdfplot(overall_shuff_acc); %%shuffled identity
h1 = cdfplot(overall_class_acc); %%pre data for post identity
h2 = cdfplot(overall_class_acc_pretest_and_prepair); %%combine all pre data
h3 = cdfplot(overall_class_acc_pretest_to_prepair); %%pretest predict prepair

h.Color = 'k';
h.LineWidth = 2;
h1.Color = 'r';
h1.LineWidth = 2;
h2.Color = 'b';
h2.LineWidth = 2;
h3.Color = [34 139 34]./255;
h3.LineWidth = 2;

legend({'Shuff', 'Pre-Post', 'AllPre', 'Day0-Day9'}, 'location', 'best')

box off
grid off

xvals = xlim;
line([xvals(1) xvals(2)], [.5 .5], 'color', [.6 .6 .6], 'linew', 1);

%%let's get stats. kstest2?
[h,p(1),ks2stat{1}] = kstest2(overall_shuff_acc, overall_class_acc);
[h,p(2),ks2stat{2}] = kstest2(overall_shuff_acc, overall_class_acc_pretest_and_prepair);
[h,p(3),ks2stat{3}] = kstest2(overall_shuff_acc, overall_class_acc_pretest_to_prepair);

saveas(gcf, fullfile(save_dir, 'CDFPlot-AccuraciesTrainOnDiffRecs.jpg'));
saveas(gcf, fullfile(save_dir, 'CDFPlot-AccuraciesTrainOnDiffRecs.svg'));