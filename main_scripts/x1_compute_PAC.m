%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of this script is to compute PAC estimates for all subjects,
% in visual timecourses from Rhodes, Rier, Singh et al. (2025)

% Author: Lyam Bailey

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define directories
restoredefaultpath
pac_config;
addpath(genpath(data_dir));
addpath(genpath(scripts_dir));
addpath(genpath(functions_dir));
addpath(fieldtrip_dir);
ft_defaults


%% Load data and define a list of subjects
all_data = load(fullfile(data_dir, "visual_gamma_data.mat"));
all_subjects_data = all_data.data;

% Define subjects
subjects = 1:102;

% Define measures
measures = {'PLV', 'ozkurt'};

% Define frequency ranges. Note that the 5-7/28-34 section of the grid cannot be used 
phase_range = [5 13];
amp_range = [28 100];

% Define time windows for the pre and post matrices (times are relative to
% stim onset)
pre_times = [-0.8 -0.1];  
post_times = [0.3 1.00];

% Seymour (2017) isolated frequencies with a 4th order butterworth filter.
% Unfortunately this can become unstable below 7 Hz. calc_MI_butter4_safe
% is a modified verion of Seymour's MI function which includes
% cfg.bpinstabilityfix = 'reduce'
filter_method = 'butter4_safe';

if strcmp(filter_method, 'butter4')
    calc_PAC = @calc_MI_orig;  % original function
elseif strcmp(filter_method, 'butter4_safe')
    calc_PAC = @calc_MI_butter4_safe;  
end


%% Loop through subjects

delete(gcp('nocreate'));  % close existing pool
parpool('local', 16);      % use 16 workers

parfor s=1:length(subjects)

    sub=subjects(s);
    
    disp('#####################')
    disp(sub)
    disp('#####################')' 

    % Point to the output folder for this subject
    subject_out_dir = fullfile(sprintf('../data/sub-%d', sub));        
    
    % Convert this subject's data to fieldtrip format
    VE_trials_ft = convert_to_fieldtrip(all_subjects_data(sub));

 
    % Compute PAC in the pre- and post-stimulus periods, using the PLV-MI
    % and Ozkurt-MI measures
    for m=1:length(measures)

        measure = measures{m};

        % Define a text string for output files, encodes frequency ranges
        % and filter method
        params_str = sprintf('%s_%d-%d_%d-%d_%s_700ms', measure, phase_range(1), phase_range(2), amp_range(1), amp_range(2), filter_method);

        % Post-stimulus period
        [matrix_post,matrix_post_surrogates] = calc_PAC(VE_trials_ft, post_times, phase_range, amp_range,'no','yes', measure);

        % Pre-stimulus period
        [matrix_pre, matrix_pre_surrogates] = calc_PAC(VE_trials_ft, pre_times, phase_range, amp_range,'no','yes', measure);

        % Save output to disks
        post_matrix_fname = fullfile(subject_out_dir, sprintf('matrix_post_%s', params_str));
        post_matrix_surrogate_fname = fullfile(subject_out_dir, sprintf('matrix_post_surrogates_%s', params_str));

        pre_matrix_fname = fullfile(subject_out_dir, sprintf('matrix_pre_%s', params_str));
        pre_matrix_surrogate_fname = fullfile(subject_out_dir, sprintf('matrix_pre_surrogates_%s', params_str));

        % A straight save() call won't work in a parfor (transparency
        % violation). Use our helper function
        args = struct();
        args.post_fname = post_matrix_fname;
        args.post_surr_fname = post_matrix_surrogate_fname;
        args.pre_fname = pre_matrix_fname;
        args.pre_surr_fname = pre_matrix_surrogate_fname;
        args.matrix_post = matrix_post;
        args.matrix_post_surrogates = matrix_post_surrogates;
        args.matrix_pre = matrix_pre;
        args.matrix_pre_surrogates = matrix_pre_surrogates;

        save_pac_results(args);


    end % measures loop

end
    
%end % subjects loop

% Helper function to save out files in the parfor loop
function save_pac_results(args)
    matrix_post = args.matrix_post;
    matrix_post_surrogates = args.matrix_post_surrogates;
    matrix_pre = args.matrix_pre;
    matrix_pre_surrogates = args.matrix_pre_surrogates;

    save(args.post_fname, 'matrix_post');
    save(args.post_surr_fname, 'matrix_post_surrogates');
    save(args.pre_fname, 'matrix_pre');
    save(args.pre_surr_fname, 'matrix_pre_surrogates');
end
