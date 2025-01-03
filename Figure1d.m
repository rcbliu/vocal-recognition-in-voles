%%%Make the stacked histogram plots for the sinFM data

clear
clc
close all

base_dir = pwd;
output_dir = 'OutputPlots';

%%sin fm feats in order within structure
sin_fm_feats = {'Duration', 'AFM', 'FFM', 'OnsetFreq', 'Slope', 'Phi'};

%%load in data
load('OrganizedVocStruc.mat');
load('VocData-tsne.mat');

addpath(genpath(fullfile(base_dir, 'Data')));
addpath(genpath(fullfile(base_dir, 'misc_code')));

%%index the data we want. Looking specifically at prepair data
num_cols = size(InfoStruc,2);
for i = 1:size(InfoStruc,1)
    phase = AllData(strcmp({AllData.AudioFname}, InfoStruc{i,1})).Phase;
    InfoStruc{i,num_cols+1} = phase;
    clear phase
end
InfoStruc = InfoStruc(find(ismember(InfoStruc(:,7), 'prepair')),:);
clear AllData

save_dir = fullfile(output_dir, 'OutputPlots\VocAnalysis');
if ~exist(save_dir)
    mkdir(save_dir);
end

%%original group colors
scatter_colors = [1,0,0;...
    0,0.800000000000000,0;...
    0,0.400000000000000,0.800000000000000;...
    0.800000000000000,0,0.800000000000000;...
    1,0.200000000000000,0.600000000000000;...
    1,0.501960784313726,0;...
    0 1 1];

%%used a website (https://mdigi.tools/lighten-color/#cc00cc) to lighten the colors. 
new_colors = [255, 153, 153;
    133, 255, 133;
    0,102,204;
    133, 194, 255;
    255, 173, 214;
    255, 204, 153;
    153, 255, 255]./255;

%%want to try to make plots where I stack all animals upon each other using
%%their individual distributions and making them semi-transparent. I don't
%%know how to describe the plot but I know how I want it to look.
%%Wish me luck.
%%step 1; get limits for each variable
sinFM_data_all = cell2mat(InfoStruc(:,5));
for i = 1:6
    min_val = nanmin(sinFM_data_all(:,i));
    max_val = nanmax(sinFM_data_all(:,i));
    bins{i,1} = linspace(min_val, max_val, 100);
    clear min_val max_val
end
clear sinFM_data_all

for i = 1:6
    fig1 = figure('color','w'); hold on
    fig2 = figure('color','w'); hold on
    yvals_additive = 0;
    
    group_data = [];
    group_labels = [];
    
    for j = 1:size(InfoStruc,1)
        data_using = InfoStruc{j,5}(:,i);
        %%now let's make our histogram
        [nelements, centers] = hist(data_using, bins{i});
        %%let's not smooth it at first
        nelements = nelements./nanmax(nelements);
        %%have to do a patch. I of course had to draw out how to do this
        %%last time I did this and I threw out those notes literally
        %%yesterday. So there's that about me.
        x_vals = cat(2,centers(1),centers,centers(end));
        
        y_vals = cat(2,yvals_additive,nelements+yvals_additive,yvals_additive);
        figure(fig1);
        patch('XData', x_vals, 'YData', y_vals, ...
            'FaceColor', new_colors(j,:), 'EdgeColor', scatter_colors(j,:));
        
        %%need to convert the data into a normal distribution
        %%log is a good choice for highly positively skewed data.
        %%May need to do separately for each variable
        %%Do need to do separately for each variable
        if j == 1
            y1 = skewness(data_using);
            if abs(y1) < 1
                disp('Normally Distributed')
                data_using_log = data_using;
                transform = 0;
            else
                fprintf('Skew = %.2f\n', y1);
                data_using_log = log10(data_using); 
                transform = 1;
            end
        else 
            if transform == 0
                data_using_log = data_using;
            else
                data_using_log = log10(data_using);
            end
        end
        
        y_vals = cat(2,yvals_additive,nelements+yvals_additive,yvals_additive);
        patch('XData', x_vals, 'YData', y_vals, ...
            'FaceColor', new_colors(j,:), 'EdgeColor', scatter_colors(j,:));
        
        yvals_additive = yvals_additive-0.5; %%Need for plotting
        group_data = cat(1,group_data,data_using_log);
        group_labels = cat(1,group_labels,ones(length(data_using),1)*j);
        clear data_using avg_val bin_diff idx_neg
        clear bin_idx height_tmp y_vals x_vals
        clear nelements centers data_using_log bins_tmp
    end
    figure(fig1);
    set(gca, 'ytick', [-3:.5:0], 'yticklabel', [7:-1:1]);
    ylabel('Rec  ID');
    xlabel(sprintf('%s', sin_fm_feats{i}));
    title_tmp = sprintf('Stacked Plot Individual sinFM feats %s', sin_fm_feats{i});

    if i == 4
        set(gca, 'xlim', [0 120000]);
    end
    
    %%let's run stats on the differences between distributions
    [p,table,stats] = anovan(group_data, group_labels);
    figure('color','w');
    [c,m,h] = multcompare(stats);
    
    %%let's get a readout of the results
    %%each row is one pairwise comparison. So let's alter our p value based
    %%on the number of comparisons
    sig_threshold = 0.05/size(c,1);
    idx_sig = find(c(:,end) <= sig_threshold);
    comparisons_sig = c(idx_sig,[1 2]);
    
    %%let's make a heatmap of sig values
    %%I think I just want the axis to go to 0.5 (even though we know
    %%significance ends well below) and color in accordingly. I'm thinking
    %%shades of red for significance, white or gray for not?
    cmap = ones(255,3);
    bins_tmp = linspace(-0.0001,0.05,254);
    bins_tmp = cat(2,bins_tmp(1)-bins_tmp(2),bins_tmp)';
    sig_data = ones(7,7)*bins_tmp(1);
    for j = 1:size(c,1)
        sig_data(c(j,1),c(j,2)) = c(j,end);
    end
    sig_bins = bins_tmp - sig_threshold;
    bin_threshold = find(sig_bins > 0);
    bin_threshold = bin_threshold(1);
    
    %%let's make the red a gradient because we can
    cmap(2:bin_threshold,1) = repmat(1,bin_threshold-1,1);
    colors_blue = linspace(0, 0.5, bin_threshold-1)';
    colors_green = colors_blue;
    cmap(2:bin_threshold,[2:3]) = cat(2,colors_blue,colors_green);
    figure('color','w'); imagesc(sig_data);
    colormap(cmap);
    caxis([bins_tmp(1) bins_tmp(end)]);
    h = colorbar;
    
    set(h, 'ytick', [0 sig_threshold 0.025 0.05]);
    h.Label.String = 'Significance';
    h.Label.Rotation = 270;
    h.Label.VerticalAlignment = "bottom";
    title_tmp2 = sprintf('Sig Heatmap VocFeatDistributions LogTransform %s', sin_fm_feats{i});
    ylabel('Pair ID1');
    xlabel('Pair ID2');
    %%let's also add the lines
    for j = 1:6
        line([0.5 7.5], [j+.5 j+.5], 'color', 'k', 'linew', 0.5);
        line([j+.5 j+.5], [0.5 7.5], 'color', 'k', 'linew', 0.5);
    end
    
    saveas(gcf, fullfile(save_dir, sprintf('%s.jpg', title_tmp2)));
    %%okay, let's instead turn the results into a table we can put into the
    %%document. Let's do the mean and standard deviation of each
    %%distribution (with the log transformed data), plus the p-value for
    %%comparison.
    %%Okay I take it back; let's do the data in a table and then do the
    %%p-values in a heatmap
    
    stat_data = zeros(7,2);
    for j = 1:length(unique(group_labels))
        data_tmp = group_data(group_labels == j);
        stat_data(j,1) = nanmean(data_tmp);
        stat_data(j,2) = nanstd(data_tmp);
    end
    
    figure(fig1);
    saveas(gcf, fullfile(save_dir, sprintf('%s.jpg', title_tmp)));
    saveas(gcf, fullfile(save_dir, sprintf('%s.svg', title_tmp)));
    sin_fm_feats{i}
    
    close all
    close all hidden
end