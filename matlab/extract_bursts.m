function bursts=extract_bursts(raw_trials, tf, times, search_freqs,...
    band_lims, fooof_thresh, sfreq, varargin)
% EXTRACT_BURSTS  Extract bursts from epoched data
%   raw_trials: raw data for each trial (trial x time)
%   tf: time-frequency decomposition for each trial (trial x freq x time)
%   times: time steps
%   search_freqs: frequency limits to search within for bursts (should be
%     wider than band_lims)
%   band_lims: keep bursts whose peak frequency falls within these limits
%   aperiodic_spectrum: aperiodic spectrum
%   sfreq: sampling rate
%   sfreq: sampling rate
%   w_size: (optional) window size to extract burst waveforms (default=0.2)
% returns: disctionary with trial, waveform, peak frequency, relative peak
%   amplitude, absolute peak amplitude, peak time, peak adjustment, FWHM in 
%   frequency, FWHM in time, and polarity for each detected burst
% Optional parameters are used as follows:
%   extract_bursts(...,'win_size',0.25)

    defaults = struct('win_size', 0.2);
    params = struct(varargin{:});
    for f = fieldnames(defaults)'
        if ~isfield(params, f{1})
            params.(f{1}) = defaults.(f{1});
        end
    end

    bursts=[];
    bursts.trial=[];
    bursts.waveform=[];
    bursts.peak_freq=[];
    bursts.peak_amp_iter=[];
    bursts.peak_amp_base=[];
    bursts.peak_time=[];
    bursts.peak_adjustment=[];
    bursts.fwhm_freq=[];
    bursts.fwhm_time=[];
    bursts.polarity=[];
    bursts.waveform_times=[];
    
    % Compute event-related signal
    erf = mean(raw_trials);

    % Iterate through trials
    for t_idx=1:size(tf,1)
        tr_tf=squeeze(tf(t_idx,:,:));

        % Regress out ERF
        lm=fitlm(erf, raw_trials(t_idx,:));
        raw_trial=lm.Residuals.Raw';

        % Extract bursts for this trial
        trial_bursts=extract_bursts_single_trial(raw_trial, tr_tf, times,...
            search_freqs, band_lims, fooof_thresh, sfreq, 'win_size',...
            params.win_size);
        
        n_trial_bursts=length(trial_bursts.peak_time);
        bursts.trial(end+1:end+n_trial_bursts)=t_idx;
        bursts.waveform(end+1:end+n_trial_bursts,:)=trial_bursts.waveform;
        bursts.peak_freq(end+1:end+n_trial_bursts)=trial_bursts.peak_freq;
        bursts.peak_amp_iter(end+1:end+n_trial_bursts)=trial_bursts.peak_amp_iter;
        bursts.peak_amp_base(end+1:end+n_trial_bursts)=trial_bursts.peak_amp_base;
        bursts.peak_time(end+1:end+n_trial_bursts)=trial_bursts.peak_time;
        bursts.peak_adjustment(end+1:end+n_trial_bursts)=trial_bursts.peak_adjustment;
        bursts.fwhm_freq(end+1:end+n_trial_bursts)=trial_bursts.fwhm_freq;
        bursts.fwhm_time(end+1:end+n_trial_bursts)=trial_bursts.fwhm_time;
        bursts.polarity(end+1:end+n_trial_bursts)=trial_bursts.polarity;
        if ~isempty(trial_bursts.waveform_times)
            bursts.waveform_times = trial_bursts.waveform_times;
        end
    end
    
end





