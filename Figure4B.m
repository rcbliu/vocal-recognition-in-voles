%%build an SVM classifier to see if we can determine a male's identity
%%based solely on his USVs.
%%Want to use USVs recorded from the pretest and prepair contexts because
%%we do not know if the USVs change based on pairing status

clear
clc
close all

%%load in data
load('OrganizedVocStruc.mat');
load('VocData-tsne.mat');
base_path = pwd;
save_dir = 'OutputPlots';
rng("default")

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));


%%randomly select an even amount of data per animal
%%Step 1: make new structure without postpair data
AllData = AllData(~strcmp({AllData.Phase}, 'postpair'));
%%then get sizes
NumVocs = [AllData.NumVocs]';

min_num_vocs = nanmin(NumVocs);

%%set up testing and training datasets
training_size = floor(min_num_vocs*0.75);
testing_size = floor(min_num_vocs*0.25);

fields = fieldnames(AllData);
fields_of_interest = [7 9:12 14 18 21];
AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);


%%Q1: Can the classifier dissociate individuals based on their pre-mating
%%USVs?
load('OrganizedVocStruc.mat');

%%randomly select an even amount of data per animal
AllData = AllData(strcmp({AllData.Phase}, 'prepair'));

AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);

tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    for j = 1:size(AllData,2)
        data_tmp = AllData(j);
        clear data_poss
        %%there is one weird USV whose duration is ridiculously short. No idea if
        %%the data got switched or what, but remove it.
        durs_tmp = AllData(j).Duration;
        idx_del = find(durs_tmp < 3);
        for k = 1:length(fields_of_interest)
            AllData(j).(fields{fields_of_interest(k)})(idx_del,:) = [];
            data_poss(:,k) = AllData(j).(fields{fields_of_interest(k)});
            data_poss(:,k) = mat2gray(data_poss(:,k));
        end
        clear randvals
        randvals = randperm(length(data_poss),min_num_vocs);
        training_data = cat(1,training_data,data_poss(randvals(1:training_size),:));
        testing_data = cat(1,testing_data,data_poss(randvals(training_size+1:training_size+testing_size),:));
        individual_ids_training(end+1:end+training_size,1) = repmat(AllData(j).Pairing,training_size,1);
        individual_ids_testing(end+1:end+testing_size,1) = repmat(AllData(j).Pairing,testing_size,1);
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
        fprintf('%2.0f%% complete', (i/1000)*100)
    end
end
toc

%%Make confusion matrix
chance_level = (100/(size(AllData,2)))/100;
figure('color','w'); hold on
[nelements,centers] = hist(prop_corr, linspace(0, 1, 50));
nelements = nelements./nansum(nelements);
h = bar(centers,nelements);
set(h,'barw',1,'facecolor','k','edgecolor','none');
yvals = ylim;
line([chance_level chance_level], [yvals(1) yvals(2)], 'linew', 2, 'color', 'r');
ylabel('Proportion of Classifiers');
xlabel('Accuracy in Determining Male');
title('Predicting Individual Identity from Pre-Mating USVs');
[h,p,ci,stats] = ttest(prop_corr, chance_level);
%%lets also add info to this plot before saving
yvals = ylim;
xvals = xlim;
%%let's assume the accuracy is greater than chance, so let's put the text
%%towards the left
text(chance_level+0.05, yvals(2)*.9, sprintf('mean = %.2f', nanmean(prop_corr)));
text(chance_level+0.05, yvals(2)*.8, sprintf('std = %.2f', nanstd(prop_corr)));
text(chance_level+0.05, yvals(2)*.7, sprintf('p = %.2f', p));
text(chance_level+0.05, yvals(2)*.6, sprintf('tstat = %.2f', stats.tstat));
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg',get(get(gca,'title'),'string'))));


confus_plot = nanmean(AllConfus,3);
cmap = flipud(gray(255));
%%make figure
figure('color','w'); imagesc(confus_plot);    
h = colorbar;
xlabel('Predicted Identity');
ylabel('True Identity');
title('Confusion Matrix - Identity Prediction Post Mating');
colormap(cmap);
caxis([0 1]);
set(h, 'xtick', [0 0.5 1]);
%%add thin black lines between boxes
yvals = ylim; xvals = xlim; %%though these should match exactly
for i = 1:length(confus_plot)
    line([xvals(1) xvals(2)], [i+.5 i+.5], 'color', 'k', 'linew', 0.25);
    line([i+.5 i+.5], [xvals(1) xvals(2)], 'color', 'k', 'linew', 0.25);
end

saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
%%add values in text
for i = 1:length(confus_plot)
    for j = 1:length(confus_plot)
        if confus_plot(i,j) > 0.6
            text(j-.2, i, sprintf('%.02f', abs(confus_plot(i,j))), 'color', 'w');
        else
            text(j-.2, i, sprintf('%.02f', abs(confus_plot(i,j))));
        end
    end
end
title('Confusion Matrix - Identity Prediction Post Mating - Add Nums');
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg',get(get(gca,'title'),'string'))));
clear cmap bins num_bins* 



%% %%Q2: Can the classifier dissociate individuals based on their pre-testing (Day 0)
%%USVs? 
%%I am finding that I can accurately dissociate the pre-testing USVs from
%%the pre-pairing USVs. Is this because there are different females?
%%Only 1 female was used for pre-testing, so if we cannot dissociate
%%identity in the pre-testing condition, maybe the female is driving our
%%results
load('OrganizedVocStruc.mat');

%%randomly select an even amount of data per animal
AllData = AllData(strcmp({AllData.Phase}, 'pretest'));

AllConfus = zeros(7,7,1000);
prop_corr = zeros(1000,1);
AllLabels = cell(1000,1);

tic
rng default
for i = 1:1000 %%make 1000 classifiers
    testing_data = [];
    training_data = [];
    individual_ids_training = [];
    individual_ids_testing = [];
    for j = 1:size(AllData,2)
        data_tmp = AllData(j);
        clear data_poss
        %%there is one weird USV whose duration is ridiculously short. No idea if
        %%the data got switched or what, but remove it.
        durs_tmp = AllData(j).Duration;
        idx_del = find(durs_tmp < 3);
        for k = 1:length(fields_of_interest)
            AllData(j).(fields{fields_of_interest(k)})(idx_del,:) = [];
            data_poss(:,k) = AllData(j).(fields{fields_of_interest(k)});
            data_poss(:,k) = mat2gray(data_poss(:,k));
        end
        clear randvals
        randvals = randperm(length(data_poss),min_num_vocs);
        training_data = cat(1,training_data,data_poss(randvals(1:training_size),:));
        testing_data = cat(1,testing_data,data_poss(randvals(training_size+1:training_size+testing_size),:));
        individual_ids_training(end+1:end+training_size,1) = repmat(AllData(j).Pairing,training_size,1);
        individual_ids_testing(end+1:end+testing_size,1) = repmat(AllData(j).Pairing,testing_size,1);
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
        fprintf('%2.0f%% complete', (i/1000)*100)
    end
end
toc

chance_level = (100/(size(AllData,2)))/100;
figure('color','w'); hold on
[nelements,centers] = hist(prop_corr, linspace(0, 1, 50));
nelements = nelements./nansum(nelements);
h = bar(centers,nelements);
set(h,'barw',1,'facecolor','k','edgecolor','none');
yvals = ylim;
line([chance_level chance_level], [yvals(1) yvals(2)], 'linew', 2, 'color', 'r');
ylabel('Proportion of Classifiers');
xlabel('Accuracy in Determining Male');
title('Predicting Individual Identity from Pretest USVs');
[h,p,ci,stats] = ttest(prop_corr, chance_level);
%%lets also add info to this plot before saving
yvals = ylim;
xvals = xlim;
%%let's assume the accuracy is greater than chance, so let's put the text
%%towards the left
text(chance_level+0.05, yvals(2)*.9, sprintf('mean = %.2f', nanmean(prop_corr)));
text(chance_level+0.05, yvals(2)*.8, sprintf('std = %.2f', nanstd(prop_corr)));
text(chance_level+0.05, yvals(2)*.7, sprintf('p = %.2f', p));
text(chance_level+0.05, yvals(2)*.6, sprintf('tstat = %.2f', stats.tstat));
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg',get(get(gca,'title'),'string'))));

confus_plot = nanmean(AllConfus,3);
bins = linspace(0, 1, 255);
cmap = flipud(gray(255));
%%make figure
figure('color','w'); imagesc(confus_plot);    
h = colorbar;
xlabel('Predicted Identity');
ylabel('True Identity');
title('Confusion Matrix - Identity Prediction Pretest');
colormap(cmap);
caxis([0 1]);
set(h, 'xtick', [0 1]);
%%add thin black lines between boxes
yvals = ylim; xvals = xlim; %%though these should match exactly
for i = 1:length(confus_plot)
    line([xvals(1) xvals(2)], [i+.5 i+.5], 'color', 'k', 'linew', 0.25);
    line([i+.5 i+.5], [xvals(1) xvals(2)], 'color', 'k', 'linew', 0.25);
end

saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg',get(get(gca,'title'),'string'))));
%%add values in text
for i = 1:length(confus_plot)
    for j = 1:length(confus_plot)
        if confus_plot(i,j) > 0.6
            text(j-.2, i, sprintf('%.02f', abs(confus_plot(i,j))), 'color', 'w');
        else
            text(j-.2, i, sprintf('%.02f', abs(confus_plot(i,j))));
        end
    end
end
title('Confusion Matrix - Identity Prediction Pretest - Add Nums');
saveas(gcf, fullfile(save_dir, sprintf('%s.jpg',get(get(gca,'title'),'string'))));
saveas(gcf, fullfile(save_dir, sprintf('%s.svg',get(get(gca,'title'),'string'))));
clear cmap bins num_bins* 