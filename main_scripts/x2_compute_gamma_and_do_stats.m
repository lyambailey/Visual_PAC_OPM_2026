%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of this script is to: (1) load in pre-computed PAC
% comodulograms (per subject), (2) compute gamma power (per subject), (3)
% get PAC clusters in children and adults, (4) perform age analyses, (5)
% make figures, (6) perform supplementary NHST analyses

% Author: Lyam Bailey
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define paths
restoredefaultpath
pac_config;
addpath(genpath(data_dir));
addpath(genpath(scripts_dir));
addpath(fieldtrip_dir);
ft_defaults

% Ensure BayesFactor is in path
addpath(genpath('../bayesFactor'));

% And colormap add-ons
addpath(genpath(fullfile('/home/lyambailey/MATLAB Add-Ons/Collections/200 colormaps/slanCM/slanCM.m')))

%% 1. Prepare Data

%% Load data

% Define list of subjects from directories in data_dir
d = dir(fullfile(data_dir, 'sub-*'));   
d = d([d.isdir]);                       
subjects = {d.name}'; 

% Sort subjects numerically
nums = cellfun(@(x) sscanf(x, 'sub-%d'), subjects);
[~, order] = sort(nums);
subjects = subjects(order);

% Drop subject 77 (repeat scan of same individual). Equivalent to 301 in
% demographics
subjects = subjects(~strcmp(subjects, 'sub-77'));

% Get count
n_subjects = length(subjects);

% Load demographic information
demographics_orig = load(fullfile(data_dir, 'demographics_alldata.mat'));
demographics = demographics_orig.demographics;

% Load original data from which PAC was computed
data_orig = load(fullfile(data_dir, "visual_gamma_data.mat")).data;

%% Analysis parameters

% Define measure of interest
measure='PLV';

% Label for figures
if strcmp(measure, 'ozkurt')
    measure_lab = 'Özkurt-MI';
elseif strcmp(measure, 'PLV')
    measure_lab = 'PLV-MI';
end

% Filter method for pulling fp and fa. 
filter_method = 'butter4_safe';  % the 'safe' version drops the filter order in cases of instability

% Define frequency ranges for PAC comodulograms
phase_range = [5 13];
amp_range = [28 100];

phase_freqs = phase_range(1):1:phase_range(2);   % 1 Hz steps
amp_freqs   = amp_range(1):2:amp_range(2);        % 2 Hz steps

% Our frequency bands are pushing reasonable ranges for PAC computation
% We'll mask out low sensitivity cells in the comodulogram(s) where
% amp_freq <= 34 Hz and phase freq >= 10 Hz
phase_min = 10;
phase_max = 13;
amp_min = 28; 
amp_max = 34;

% Define time windows for pre/post stimulus
pre_times = [-0.8 -0.1];  
post_times = [0.3 1.00];

% Define the range of broadband gamma
gamma_range = [30 80];
gamma_freqs = gamma_range(1):2:gamma_range(2);  % 2 Hz steps


%% Load PAC estimates for all subjects and compute gamma power

% Define and empty struct to store PAC and gamma power estimates (pre vs
% post) and demographic info for all subjects
clear results
results(length(subjects)) = struct('subjectID', [], 'age', [], 'sex', [], ...
    'pac_pre', [], ...  % pre-stimulus pac estimates
    'pac_post', [], ... % post-stimulus pac estimates
    'pac_diff', [], ... % post minus pe-stimulus pac
    'gamma_rel_change', [], ...  % relative gamma [30-80 Hz] change between pre and post, at every 2 Hz frequency step 
    'loaded', false);  % indicates whether subject data failed to load

% Loop through subjects
for s = 1:length(subjects)
    subject = subjects{s};

    fprintf('Iter %d: subject=%s, subject_idx=%d\n', s, subject, sscanf(subject, 'sub-%d'))

    % Subjects are labeled according to their index in all_data (original
    % data from which PAC was computed). Entries in demographics are in the
    % same order. So, get the index from the subject ID, and then use that
    % to pull their demographic information
    subject_idx = sscanf(subject, 'sub-%d');

    % Pull this subject's demographics data
    subject_age = demographics.age(subject_idx);
    subject_sex = demographics.sex(subject_idx);

    % Compute broadband gamma change...

    % Pull the VE for this subject. As with demographics, index with
    % subject_idx
    VE_trials_ft = convert_to_fieldtrip(data_orig(subject_idx));

    % Get gamma power in each time window    
    [fa_post, ~, mean_gammas_post] = get_gamma_power(VE_trials_ft, gamma_range, post_times);
    [fa_pre, ~, mean_gammas_pre] = get_gamma_power(VE_trials_ft, gamma_range, pre_times);

    % Compute relative change: (post-pre) / pre
    gamma_rel_change = (mean_gammas_post - mean_gammas_pre) ./ mean_gammas_pre;  % one value per frequency step
    
    % Load pac estimates
    try
        % Load pre- and post-stimulus PAC matrices (surrogate corrected) for this subject
        params_str = sprintf('%s_%d-%d_%d-%d_%s_700ms', measure, phase_range(1), phase_range(2), amp_range(1), amp_range(2), filter_method);

        post = load(fullfile(data_dir, subject, sprintf('matrix_post_surrogates_%s.mat', params_str)));
        pre = load(fullfile(data_dir, subject, sprintf('matrix_pre_surrogates_%s.mat', params_str)));

        % Compute the difference between post and pre
        pac_diff = post.matrix_post_surrogates - pre.matrix_pre_surrogates;

        % OR use relative change (as we did with gamma)
        pac_rel_change =  (post.matrix_post_surrogates - pre.matrix_pre_surrogates) ./ pre.matrix_pre_surrogates;
        
        % Store demographic info and results for this subject.
        results(s).subjectID = subject;
        results(s).age = subject_age;
        results(s).sex = subject_sex;
        results(s).pac_pre = pre.matrix_pre_surrogates;
        results(s).pac_post = post.matrix_post_surrogates;
        results(s).pac_diff =  pac_diff;
        results(s).gamma_rel_change = gamma_rel_change;  % note that this stores gamma change for all freqs. 

        % Mark this subject as having successfully loaded data
        results(s).loaded = true;

    catch ME
        fprintf('Skipping %s: %s\n', subject, ME.message);
        results(s).loaded = false;
    
    end % try-catch loop
end % subjects loop

% Drop empty slices (i.e., where we failed to load data for a given
% subject
n_dropped = n_subjects - sum([results.loaded]);
fprintf('dropping %d empty slices', n_dropped);
results = results([results.loaded]);

% Ensure that each subject's age in results still matches its entry in
% demographics
idx = arrayfun(@(r) sscanf(r.subjectID, 'sub-%d'), results);
assert(isequal([results.age]', demographics.age(idx)), ...
    'Age mismatch between results and demographics');

% Save out results
save(fullfile(results_dir, sprintf('%s_pac_results_all_subjects.mat', measure)), 'results');

%% Optionally load pre-computed results
results = load(fullfile(results_dir, sprintf('%s_pac_results_all_subjects.mat', measure))).results;

%% Extract age, PAC, and gamma arrays from results
ages = [results.age]';

% Concatenate PAC matrices across subjects
pac_change = cat(3, results.pac_diff);  % amp_freqs x phase_freqs x subjects array

% Get mean broadband gamma per subject
gamma_array = cat(1, results.gamma_rel_change);   % subjects x frequencies array
gamma_broad = mean(gamma_array(:, gamma_freqs >= 30 & gamma_freqs <= 80), 2);  % one value per subject

% Apply rank transforms to all the above
ages_rank = tiedrank(ages);
pac_diffs_rank = tiedrank(pac_change);
gamma_broad_rank = tiedrank(gamma_broad);

%% 2. Get clusters for pre vs post-stim PAC
% We'll test for the presence of PAC independently in each age group

% Define age groups and subset data for each
is_child =  ages < 18;   
is_adult = ages >= 18;

results_child = results(is_child);
results_adult = results(is_adult);
matrices_pre_child = cat(3, results_child.pac_pre);
matrices_post_child = cat(3, results_child.pac_post);
matrices_pre_adult = cat(3, results_adult.pac_pre);
matrices_post_adult = cat(3, results_adult.pac_post);

% Compute grand average per age group
grand_child = mean(pac_change(:, :, is_child), 3);
grand_adult = mean(pac_change(:, :, is_adult), 3);

% Create empty arrays to store BFs from children and adults
pac_dims = size(results(1).pac_pre); % pre, post, and diff should all be equivalent size
bf_map_child = zeros(pac_dims);
bf_map_adult = zeros(pac_dims);

% Loop through matrix rows (amplitude freqs)
for j = 1:pac_dims(1)
    % Loop through matrix columns (phase freqs)
    for k = 1:pac_dims(2)

        % Get pre and post values for all subjects at the jth/kth position,
        % independently per group
        pre_values_child = squeeze(matrices_pre_child(j, k, :));   
        post_values_child = squeeze(matrices_post_child(j, k, :));   

        pre_values_adult = squeeze(matrices_pre_adult(j, k, :));   
        post_values_adult = squeeze(matrices_post_adult(j, k, :)); 

        % Do a paired-samples Bayes t tests between pre and post
        bf10_prepost_child = bf.ttest(post_values_child, pre_values_child, 'Tail', 'right');
        bf10_prepost_adult = bf.ttest(post_values_adult, pre_values_adult, 'Tail', 'right');
        
        % Store BFs
        bf_map_child(j,k) = bf10_prepost_child;
        bf_map_adult(j,k) = bf10_prepost_adult;

    end % k loop
end % j loop

% Create a logical mask from the top 5% of values across both (child/adult)
% BF maps.

% Threshold values in each map
threshold = quantile([bf_map_child(:); bf_map_adult(:)], 0.95);
assert(threshold > 1.0, "Warning, the 5% threshold is < 1.0");  % threshold should exceed 1.0

mask_top5_child = bf_map_child >= threshold;
mask_top5_adult = bf_map_adult >= threshold;

% Find clusters in each map
CC_adult = bwconncomp(mask_top5_adult);
CC_child = bwconncomp(mask_top5_child);

% Get cluster sizes, and find the largest cluster in each map
sizes_adult = cellfun(@numel, CC_adult.PixelIdxList);
sizes_child = cellfun(@numel, CC_child.PixelIdxList);
[max_adult, idx_adult] = max(sizes_adult);
[max_child, idx_child] = max(sizes_child);

% Take whichever cluster is largest overall
if max_adult >= max_child
    mask_largest = false(size(mask_top5_adult));
    mask_largest(CC_adult.PixelIdxList{idx_adult}) = true;
else
    mask_largest = false(size(mask_top5_child));
    mask_largest(CC_child.PixelIdxList{idx_child}) = true;
end

% Get linear indices of the chosen cluster
if max_adult >= max_child
    cluster_lin = CC_adult.PixelIdxList{idx_adult};
else
    cluster_lin = CC_child.PixelIdxList{idx_child};
end

% Convert to row (amp) and column (phase) indices
[rows, cols] = ind2sub(size(mask_largest), cluster_lin);

% Map to actual frequencies
cont_min_phase = phase_freqs(min(cols));
cont_max_phase = phase_freqs(max(cols));
cont_min_amp   = amp_freqs(min(rows));
cont_max_amp   = amp_freqs(max(rows));

% Print out cluster ranges
fprintf('Cluster spans phase %g–%g Hz, amplitude %g–%g Hz\n', ...
    cont_min_phase, cont_max_phase, cont_min_amp, cont_max_amp);

% For each subject, compute the mean PAC difference within the largest
% cluster. 
cluster_pac_change = zeros(n_subjects, 1);

for s = 1:n_subjects
    % Get this subject's difference matrix
    slice = pac_change(:, :, s);

    % Get their mean value within the mask
    cluster_pac_change(s) = mean(slice(mask_largest));

end

%% Main Stats

% Perform Bayesian covariate testing. To get a BF for the "main effect" of
% some variable (age or gamma), we need to compute two models - one full
% model (including that variable) and one restricted (including only the
% covariates), and then get the BF for the model comparison.

% Create data table for the regressions below. The different measures are
% on very different scales, which the JZS prior struggles with. We'll
% zscore all measures
tbl = table(cluster_pac_change, ages, gamma_broad, ...
    'VariableNames', {'pac','age','gamma'});

tbl_z = varfun(@zscore, tbl);

% Fit models
model_full = bf.anova(tbl_z, 'zscore_pac ~ zscore_age + zscore_gamma', 'Verbose', 0); % includes age + gamma
model_no_age = bf.anova(tbl_z, 'zscore_pac ~ zscore_gamma', 'Verbose', 0); % includes gamma only
model_no_gamma = bf.anova(tbl_z, 'zscore_pac ~ zscore_age', 'Verbose', 0); % includes gamma only

% Compare models to approximate "main effect of age"
bf_age = model_full / model_no_age;
bf_gamma = model_full / model_no_gamma;
disp(sprintf('Age effect BF: %f', bf_age));
disp(sprintf('Gamma effect BF: %f', bf_gamma));


% Non-bayesian multiple regression on ranks
cluster_pac_change_rank = tiedrank(cluster_pac_change);
tbl_rnk = table(cluster_pac_change_rank, ages_rank, gamma_broad_rank, ...
    'VariableNames', {'pac','age','gamma'});

tbl_rnk_z = varfun(@zscore, tbl_rnk);

lm = fitlm(tbl_rnk, 'pac ~ age + gamma');
disp(lm.Coefficients);


%% Plot pre/post GAVG for each age group with stat maps

% Define shared color scale for GAVG plots and stat maps
cmax = max(abs([grand_child(:); grand_adult(:)]))*0.9;
logbmax = 1.7;  % equal to BF ~= 50 
bticks = 0:0.25:ceil(logbmax * 10)/10; % rounds logbmax up to one decimal place

% Set up figure
figure;
nrow = 2;
ncol = 2;
t = tiledlayout(nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top row: child and adult GAVG's

% Child GAVG
ax1 = nexttile(1);
pcolor(phase_freqs, amp_freqs, grand_child);
shading interp;
axis xy; colormap; clim([-cmax cmax]);
%xlabel('Phase frequency (Hz)'); 
ylabel('Amplitude frequency (Hz)');
title('Children');
cb = colorbar(ax1);
cb.Visible = 'off';

% Overlay white rectangle for low sensitivity zone
hold on;
rectangle('Position', [phase_min amp_min phase_max-phase_min amp_max-amp_min], ...
    'FaceColor', '#D3D3D3', 'EdgeColor', 'none');
hold off;

% Adult GAVG
ax2 =  nexttile(2);
pcolor(phase_freqs,amp_freqs, grand_adult);
shading interp;
axis xy; colormap; clim([-cmax cmax]);
cb = colorbar(ax2);
cb.Label.String = 'ΔPAC'; 
cb.Ticks = [-cmax, cmax];

cb.TickLabels = {'-', '+'};
cb.FontSize=12;
cb.Label.FontSize = 12;
%xlabel('Phase frequency (Hz)'); 
%ylabel('Amplitude frequency (Hz)');
title('Adults');

hold on;
rectangle('Position', [phase_min amp_min phase_max-phase_min amp_max-amp_min], ...
    'FaceColor', '#D3D3D3', 'EdgeColor', 'none');

% Bottom row: child and adult stat maps
ax3 =  nexttile(3);
pcolor(phase_freqs, amp_freqs, log10(bf_map_child));
axis xy;
clim([0 logbmax]);     
cb = colorbar(ax3);
cb.Visible = 'off';  
shading interp;
xlabel('Phase frequency (Hz)');
ylabel('Amplitude frequency (Hz)');

% Overlay white rectangle for low sensitivity zone
hold on;
rectangle('Position', [phase_min amp_min phase_max-phase_min amp_max-amp_min], ...
    'FaceColor', '#D3D3D3', 'EdgeColor', 'none');
hold off;

ax4 =  nexttile(4);
pcolor(phase_freqs, amp_freqs, log10(bf_map_adult));
axis xy;
clim([0 logbmax]);     
cb = colorbar;
cb.Ticks = bticks;
cb.TickLabels = arrayfun(@(x) sprintf('%.d', round(10^x)), bticks, 'UniformOutput', false); 
cb.Label.String = 'BF_{Stimulus}';
cb.Label.FontSize = 12;
shading interp;
xlabel('Phase frequency (Hz)');
%ylabel('Amplitude frequency (Hz)');

% Plot cluster on adult map
hold on;
contour(phase_freqs, amp_freqs, double(mask_largest), 1, 'k-', 'LineWidth', 2); % show the top BF cluster
hold off;

% Overlay white rectangle for low sensitivity zone
hold on;
rectangle('Position', [phase_min amp_min phase_max-phase_min amp_max-amp_min], ...
    'FaceColor', '#D3D3D3', 'EdgeColor', 'none');
hold off;

% Set colormaps
colormap(ax1, slanCM('vik'));
colormap(ax2, slanCM('vik'));
colormap(ax3, slanCM('amp'));
colormap(ax4, slanCM('amp'));

% Add a/b/c/d to each panel
axes_list = [ax1, ax2, ax3, ax4];
labels    = {'A', 'B', 'C', 'D'};
for i = 1:numel(axes_list)
    text(axes_list(i), -0.15, 1.05, labels{i}, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 16, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

% save figure
exportgraphics(gcf, fullfile(results_dir, sprintf('%s_pac_adults_vs_kids_gavg.png', measure)), 'Resolution', 300);


%% Plot the PAC ~ age relationship

figure;
scatter(ages(is_child), cluster_pac_change(is_child));
%ylim([-4,3]);
lsline;

% Add stats to figure. Note that in the lm coefficients table, variables
% are indexed by row; stats by column. 
% Row 1 = Intercept, 2 = Age, 3 = Gamma.
b = lm.Coefficients(2,1).Estimate;  % B estimate for age
p = lm.Coefficients(2,4).pValue; % p value for age

text(0.65, 0.25, sprintf('BF_{Age} = %.2f\n\\beta_{rank} = %.2f\n\\itp\\rm_{rank} = %.2f', ...
    bf_age, b, p), ...
    'Units', 'normalized', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'FontSize', 10, 'Interpreter', 'tex', ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');
ylabel('Mean Δ PAC')
xlabel('Age')
title(sprintf('Age-PAC relationship within largest cluster \n [Phase: %g-%g Hz, Amp: %g-%g Hz]', ...
    cont_min_phase, cont_max_phase, cont_min_amp, cont_max_amp));


% Save figure
set(gcf, 'Units', 'inches', 'Position', [1 1 5 4]);   % 5×5 inch figure
exportgraphics(gcf, fullfile(results_dir, sprintf('%s_pac_age_scatter.png', measure)), 'Resolution', 300);


%% Standard stats, for supplemental analysis
kids = {results_child.subjectID};
adults = {results_adult.subjectID};

params_str = sprintf('%s_%d-%d_%d-%d_%s_700ms', measure, phase_range(1), phase_range(2), amp_range(1), amp_range(2), filter_method);

[stat_child] = get_PAC_stats_righttailed(sprintf('matrix_post_surrogates_%s.mat',params_str),...
    sprintf('matrix_pre_surrogates_%s', params_str),[5 13],[28 100],kids,'../data/',1);
make_smoothed_comodulograms(stat_child, [5 13], [28 100]);
colormap(slanCM('vik'));
 
[stat_adult] = get_PAC_stats_righttailed(sprintf('matrix_post_surrogates_%s.mat',params_str),...
    sprintf('matrix_pre_surrogates_%s', params_str),[5 13],[28 100],adults,'../data/',1);
make_smoothed_comodulograms_fixedclim(stat_adult, [5 13], [28 100], [-3.25, 3.25]);
colormap(slanCM('vik'));

% Add title
if strcmp(measure, 'ozkurt')
    title_lab = 'Özkurt-MI';
elseif strcmp(measure, 'PLV')
    title_lab = 'PLV-MI';
end

title(title_lab, 'FontSize', 20, 'FontWeight', 'bold');

set(gcf, 'Units', 'inches', 'Position', [0.5 0.5 10 10]);   % 10x10 inch figure
exportgraphics(gcf, fullfile(results_dir, sprintf('%s_pac_cluster_stats.png', measure)), 'Resolution', 300);

% Get linear indices of the significant (adult) cluster
cluster_lin = find(squeeze(stat_adult.mask));
[rows, cols] = ind2sub(size(squeeze(stat_adult.mask)), cluster_lin);

cont_min_phase = phase_freqs(min(cols));
cont_max_phase = phase_freqs(max(cols));
cont_min_amp   = amp_freqs(min(rows));
cont_max_amp   = amp_freqs(max(rows));

fprintf('Cluster spans phase %g–%g Hz, amplitude %g–%g Hz\n', ...
    cont_min_phase, cont_max_phase, cont_min_amp, cont_max_amp);

% Constrain subjects' comodulogram to the significant cluster, compute mean
% within the cluster and correlate with age (similar to above)
mean_pac_diff_clusterstats = zeros(n_subjects, 1);

for s = 1:n_subjects
    % Get this subject's difference matrix
    slice = pac_change(:, :, s);

    % Get their mean value within the mask
    mean_pac_diff_clusterstats(s) = mean(slice(stat_adult.mask));
end

[rho, p_val] = corr(ages, mean_pac_diff_clusterstats, 'Type','Spearman');

% Apply Bonferroni correction
p_val_bonf = p_val * length(measures);

% Show scatter
figure;
scatter(ages, zscore(mean_pac_diff_clusterstats));
ylim([-4,3]);
lsline;

% Add stats to figure
text(0.78, 0.15, sprintf('\\itr\\rm_{S} = %.3f\n\\itp\\rm = %.3f', rho, p_val_bonf), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 11, 'BackgroundColor', 'w', 'EdgeColor', 'k');
ylabel('z-scored PAC change (post - pre)')
xlabel('Age')
title(sprintf('%s: Age-PAC correlation within significant cluster \n [Phase: %g-%g Hz, Amp: %g-%g Hz]', ...
    measure_lab, cont_min_phase, cont_max_phase, cont_min_amp, cont_max_amp));


% Save figure
set(gcf, 'Units', 'inches', 'Position', [1 1 5 5]);   % 5×5 inch figure
exportgraphics(gcf, fullfile(results_dir, sprintf('%s_pac_age_scatter_from_cluster_stats.png', measure)), 'Resolution', 300);