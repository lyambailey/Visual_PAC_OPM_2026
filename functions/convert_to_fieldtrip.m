function VE_ft = convert_to_fieldtrip(d)
% Helper function to convert VE structs into Fieldtrip-compatible struct, suitable for calc_MI.

fs          = d.samp_frequency;          % 1200 Hz
n_samples   = size(d.VE_trials, 1);      % 3600
n_trials    = d.N_trials;                % 51
trig_offset = d.trig_offset;             % -1 s

t = (0:n_samples-1) / fs + trig_offset; % -1 to +2 s

VE_ft.label      = {'VE'};
VE_ft.fsample    = fs;
VE_ft.trialinfo  = (1:n_trials)';

for k = 1:n_trials
    VE_ft.trial{1,k} = d.VE_trials(:, k)';  % [1 x 3600]
    VE_ft.time{1,k}  = t;                    % [1 x 3600]
end

% sampleinfo: [start_sample, end_sample] for each trial
VE_ft.sampleinfo = [(0:n_trials-1)' * n_samples + 1, ...
                    (1:n_trials)'   * n_samples];
end