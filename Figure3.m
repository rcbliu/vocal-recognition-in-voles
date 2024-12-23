%%Want to use the framework I have for plotting the tsne and plotting pre
%%vs post-mating across all adults and instead break it down by individual.

clear
clc
close all

make_plots = 0; %%make the plots? yes (1) or no (0).

scatter_colors = [1,0,0;0,0.800000000000000,0;
    0,0.400000000000000,0.800000000000000;
    0.800000000000000,0,0.800000000000000;
    1,0.200000000000000,0.600000000000000;
    1,0.501960784313726,0;
    0.627450980392157,0.627450980392157,0.627450980392157];

base_path = pwd;
save_dir_figures = 'OutputPlots';

addpath(genpath(fullfile(base_path, 'Data')));
addpath(genpath(fullfile(base_path, 'misc_code')));

load('VocData-tsne');

prepair_names = {'ch1T0000626', ...
    'ch1T0000627', 'ch1T0000628', 'ch1T0000629', 'ch1T0000642', ...
    'ch1T0000643', 'ch1T0000644'};

postpair_names = {'ch1T0000638', 'ch1T0000639', 'ch1T0000640', ...
    'ch1T0000641', 'ch1T0000677', 'ch1T0000678', 'ch1T0000679'};

var_labels = {'Duration', 'Amp', 'Freq', 'OnsetFreq', 'Slope', 'Phi'};

num_files = length(prepair_names)+length(postpair_names);

%%Organize data into precohab and postcohab
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

%%get tsne data
all_y_data = cat(1,y_data_prepair,y_data_postpair);
labels_y(1:length(data_prepair),1) = 1;
labels_y(end+1:end+length(data_postpair),1) = 2;

%%run density function on tsne data
min_val = nanmin(nanmin(all_y_data));
max_val = nanmax(nanmax(all_y_data));
fig_all = figure('color','w');
[xx_postpair,density_postpair] = findPointDensity(all_y_data, ...
    3,1001,[min_val max_val]);
if make_plots == 1
    imagesc(density_postpair); colormap(jet(255));
    xvals = xlim;
    set(gca, 'xtick', linspace(xvals(1), xvals(2), 11), 'ytick', linspace(xvals(1), xvals(2),11));
    set(gca, 'xticklabel', [-50:10:50], 'yticklabel', [-50:10:50]);
    set(gca, 'ydir', 'normal');
    title_name = 'tSNE Density Gaussian AllVocs';
    title(title_name);
    axis square
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap AllVocs.jpg'));
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap AllVocs.svg'));
end

all_vals = cell2mat(InfoStruc(:,6));

%%let's make comparisons for all pre vs all post
data_prepair_tmp = cell2mat(cell_data_prepair(:,2));
data_postpair_tmp = cell2mat(cell_data_postpair(:,2));
[xx_prepair,density_prepair] = findPointDensity(data_prepair_tmp, ...
    3,1001,[min_val max_val]);
[xx_postpair,density_postpair] = findPointDensity(data_postpair_tmp, ...
    3,1001,[min_val max_val]);
dist_pre_post = nansum(JSDiv(density_prepair, density_postpair))

%%now make plots organized by individual
for i = 1:length(cell_data_prepair)
    data_pre_tmp = cell_data_prepair{i,2};
    data_post_tmp = cell_data_postpair{i,2};
    
    min_val = nanmin(nanmin(cat(1,data_pre_tmp,data_post_tmp)));
    max_val = nanmax(nanmax(cat(1,data_pre_tmp,data_post_tmp)));
    
    [xx_prepair,density_prepair] = findPointDensity(data_pre_tmp, ...
        3,1001,[min_val max_val]);
    [xx_postpair,density_postpair] = findPointDensity(data_post_tmp, ...
        3,1001,[min_val max_val]);
    
    if make_plots == 1
        fig_pre = figure('color','w');
        imagesc(density_prepair); colormap(jet(255));
        xvals = xlim;
        set(gca, 'xtick', linspace(xvals(1), xvals(2), 11), 'ytick', linspace(xvals(1), xvals(2),11));
        set(gca, 'xticklabel', [-50:10:50], 'yticklabel', [-50:10:50]);
        set(gca, 'ydir', 'normal');
        title_name = sprintf('tSNE Density Gaussian Prepair %i', i);
        title(title_name);
        axis square
        
        saveas(gcf, fullfile(save_dir_figures, sprintf('Tsne Heatmap Pre - Pair %i.jpg', i)));
        saveas(gcf, fullfile(save_dir_figures, sprintf('Tsne Heatmap Pre - Pair %i.svg', i)));
        
        fig_post = figure('color','w');
        
        imagesc(density_postpair); colormap(jet(255));
        xvals = xlim;
        set(gca, 'xtick', linspace(xvals(1), xvals(2), 11), 'ytick', linspace(xvals(1), xvals(2),11));
        set(gca, 'xticklabel', [-50:10:50], 'yticklabel', [-50:10:50]);
        set(gca, 'ydir', 'normal');
        title_name = sprintf('tSNE Density Gaussian Postpair %i', i);
        title(title_name);
        axis square
        saveas(gcf, fullfile(save_dir_figures, sprintf('Tsne Heatmap Post - Pair %i.jpg', i)));
        saveas(gcf, fullfile(save_dir_figures, sprintf('Tsne Heatmap Post - Pair %i.svg', i)));
    end
    
    all_density{i,1} = density_prepair;
    all_density{i,2} = density_postpair;
    
    clear data_*_tmp min_val max_val
    clear density_prepair density_postpair xx_*
end

%%now let's get the distances between every possible pair of recordings
dists_pre_post = NaN(length(all_density),1);
dists_individs = NaN(length(all_density), length(all_density),2);
for i = 1:length(all_density)
    dists_pre_post(i,1) = nansum(JSDiv(all_density{i,1},all_density{i,2}))/100;
    for j = i+1:length(all_density)
        dists_individs(i,j,1) = nansum(JSDiv(all_density{i,1}, all_density{j,1}))/100;
        dists_individs(i,j,2) = nansum(JSDiv(all_density{i,2}, all_density{j,2}))/100;
    end
end

max_dist = ceil(max(max(max(dists_individs))));

%%dist_vals_all(:,1) = pre-to-post distance within an individual
%%dist_vals_all(:,2:end) = distances between pairwise individual
%%comparisons
for i = 1:length(dists_individs)
    count = 2;
    for j = 1:length(dists_individs)
        if ~isnan(dists_individs(i,j,1))
            dist_vals_all(i,count) = dists_individs(i,j,1);
            dist_vals_all(i,count+1) = dists_individs(i,j,2);
            count = count+2;
        elseif ~isnan(dists_individs(j,i,1))
            dist_vals_all(i,count) = dists_individs(j,i,1);
            dist_vals_all(i,count+1) = dists_individs(j,i,2);
            count = count+2;
        end
    end
end
dist_vals_all(:,1) = dists_pre_post;

if make_plots == 1    
    %%run stats
    %%find the distribution of difference values and see if it differs from
    %%zero
    dist_vals_tmp = [];
    for i = 1:size(dist_vals_all,1)
        dist_vals_tmp = cat(2,dist_vals_tmp,dist_vals_all(i,2:end)-dist_vals_all(i,1));
    end
    
    z_test = (0 - nanmean(dist_vals_tmp))./nanstd(dist_vals_tmp);
    
    %%run a ttest for within versus between values
    btwn_vals = reshape(dist_vals_all(:,2:end),1,size(dist_vals_all,1)*(size(dist_vals_all,2)-1));
    [h,p] = ttest2(dist_vals_all(:,1)',btwn_vals)

    %%also make a plot
    figure('color','w'); hold on
    h = cdfplot(btwn_vals);
    set(h, 'color', 'k', 'linew', 2);
    h = cdfplot(dists_pre_post);
    set(h, 'color', 'r', 'linew', 2);
    saveas(gcf, fullfile(save_dir_figures, 'CDFPlot-DivergenceWithinVsBtwn.tiff'));
    saveas(gcf, fullfile(save_dir_figures, 'CDFPlot-DivergenceWithinVsBtwn.svg'));
    
    %%Run a ttest to compare the difference between all between-animal within-context and
    %%their corresponding within-animal between-context values to zero
    [h,p] = ttest(dist_vals_tmp)

    %%Find whether each value is above or below the
    %%'within-vole' value and run a sign rank
    foo = dist_vals_all(:,2:end)-dist_vals_all(:,1);
    for i = 1:size(foo,1)
        num_pos(i,1) = length(find(foo(i,:) > 0));
    end
    prop_pos = num_pos./length(foo);
    vals = cat(2,num_pos,length(foo)-num_pos);
    vals(:,3) = length(foo)/2;
    
    foo = cat(1,sum(num_pos),sum(vals(:,2)));
    foo(:,2) = sum(foo)/2;
    foo(:,3) = (foo(:,2)-foo(:,1)).^2./foo(:,2);
    chi_sq = sum(foo(:,3))
    p=1-chi2cdf(chi_sq,1)

    %%set up data for plot
    normalized_dist = bsxfun(@minus,dist_vals_all(:,2:end),dist_vals_all(:,1));
    normalized_dist = reshape(normalized_dist,size(normalized_dist,1)*size(normalized_dist,2),1);
    figure('color','w'); h = cdfplot(normalized_dist); hold on
    set(h, 'color', 'k', 'linew', 2);
    yvals = ylim;
    line([0 0], [yvals(1) yvals(2)], 'color', 'r', 'linew', 2, 'linestyle', '--');
    iqr = prctile(normalized_dist, [25 50 75]);
    line([iqr(1) iqr(3)], [0.5 0.5], 'color', [.7 .7 .7], 'linew', 2);
    line([iqr(1) iqr(1)], [0.45 0.55], 'color', [.7 .7 .7], 'linew', 2);
    line([iqr(3) iqr(3)], [0.45 0.55], 'color', [.7 .7 .7], 'linew', 2);
    scatter(iqr(2), 0.5, 60, 'b', 'filled');
    saveas(gcf, fullfile(save_dir_figures, 'AxMinusWithinDivergence.tiff'));
    saveas(gcf, fullfile(save_dir_figures, 'AxMinusWithinDivergence.svg'));
    
    %%now do prepair and postpair as separate plots
    all_data_pre = [];
    all_data_post = [];
    for i = 1:size(cell_data_prepair,1)
        all_data_pre = cat(1,all_data_pre,cell_data_prepair{i,2});
        all_data_post = cat(1,all_data_post,cell_data_postpair{i,2});
    end
    min_val = nanmin(nanmin(all_data_pre));
    max_val = nanmax(nanmax(all_data_pre));
    
    fig_pre = figure('color','w');
    [xx_prepair,density_prepair] = findPointDensity(all_data_pre, ...
        3,1001,[min_val max_val]);
    imagesc(density_prepair); colormap(jet(255));
    xvals = xlim;
    set(gca, 'xtick', linspace(xvals(1), xvals(2), 11), 'ytick', linspace(xvals(1), xvals(2),11));
    set(gca, 'xticklabel', [-50:10:50], 'yticklabel', [-50:10:50]);
    set(gca, 'ydir', 'normal');
    title_name = 'Tsne Heatmap All Adult Pre Data';
    title(title_name);
    axis square
    
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap All Adult Pre Data.jpg'));
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap All Adult Pre Data.svg'));
    
    min_val = nanmin(nanmin(all_data_post));
    max_val = nanmax(nanmax(all_data_post));
    
    fig_pre = figure('color','w');
    [xx_prepair,density_prepair] = findPointDensity(all_data_post, ...
        3,1001,[min_val max_val]);
    imagesc(density_prepair); colormap(jet(255));
    xvals = xlim;
    set(gca, 'xtick', linspace(xvals(1), xvals(2), 11), 'ytick', linspace(xvals(1), xvals(2),11));
    set(gca, 'xticklabel', [-50:10:50], 'yticklabel', [-50:10:50]);
    set(gca, 'ydir', 'normal');
    title_name = 'Tsne Heatmap All Adult Post Data';
    title(title_name);
    axis square
    
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap All Adult Post Data.jpg'));
    saveas(gcf, fullfile(save_dir_figures, 'Tsne Heatmap All Adult Post Data.svg'));
end
