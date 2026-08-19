%% Function for getting pre/post gamma change
function [fa, mean_gammas_trials, mean_gammas] = get_gamma_power(VE, freq_range, toi)

% mean_gammas_trials is a trials x fa array of amplitude estimates, where
% fa is a list of center frequencies defined by amp_range. mean_gammas is a
% 1 x fa array, averaged over trials

fa = freq_range(1):2:freq_range(2);

% Pre-allocate empty array to store mean gamma amplitude over epochs
mean_gammas_trials = nan(numel(VE.trial), numel(fa));

for i = 1:numel(fa)

    % Apply two-pass, fourth-order butterworth filter (center freq +/- 2 Hz)
    cfg = [];
    cfg.showcallinfo = 'no';
    cfg.bpfilttype = 'but';
    cfg.bpfiltord = 4;
    cfg.bpfiltdir = 'twopass';
    cfg.hilbert = 'abs';
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [fa(i)-2 fa(i)+2];

    VE_gamma = ft_preprocessing(cfg, VE);

    % Split into pre- and post-stimulus
    cfg = [];
    cfg.showcallinfo = 'no';
    cfg.toilim = toi;
    VE_gamma_toi = ft_redefinetrial(cfg, VE_gamma);

    % Average over timepoints and trials
    mean_gammas_trials(:, i) = cellfun(@(x) mean(x(1,:)), VE_gamma_toi.trial(:));

end % amplitude loop

% Average over trials
mean_gammas = mean(mean_gammas_trials, 1);

end % function