function bursts=extract_bursts_single_trial(raw_trial, tf, times,...
    search_freqs, band_lims, aperiodic_spectrum, sfreq, varargin)
% EXTRACT_BURSTS_SINGLE_TRIAL  Extract bursts from epoched data
%   raw_trials: raw data for trial (time)
%   tf: time-frequency decomposition for trial (freq x time)
%   times: time steps
%   search_freqs: frequency limits to search within for bursts (should be
%     wider than band_lims)
%   band_lims: keep bursts whose peak frequency falls within these limits
%   aperiodic_spectrum: aperiodic spectrum
%   sfreq: sampling rate
%   sfreq: sampling rate
%   w_size: (optional) window size to extract burst waveforms (default=0.2)
% returns: disctionary with waveform, peak frequency, relative peak
%   amplitude, absolute peak amplitude, peak time, peak adjustment, FWHM in 
%   frequency, FWHM in time, and polarity for each detected burst
% Optional parameters are used as follows:
%   extract_bursts(...,'win_size',0.25)

    defaults = struct('win_size', .2);
    params = struct(varargin{:});
    for f = fieldnames(defaults)'
        if ~isfield(params, f{1})
            params.(f{1}) = defaults.(f{1});
        end
    end

    bursts=[];
    bursts.waveform=[];
    bursts.peak_freq=[];
    bursts.peak_amp_iter=[];
    bursts.peak_amp_base=[];
    bursts.peak_time=[];
    bursts.peak_adjustment=[];
    bursts.fwhm_freq=[];
    bursts.fwhm_time=[];
    bursts.polarity=[];
    bursts.waveform_times = [];
    
    % Grid for computing 2D Gaussians
    [x_idx, y_idx] = meshgrid(1:length(times), 1:length(search_freqs));

    % Window size in points
    wlen = round(params.win_size * sfreq);
    half_wlen = round(wlen * .5);

    % Subtract 1/f threshold
    trial_tf = tf - repmat(aperiodic_spectrum,1,size(tf,2));
    trial_tf(trial_tf < 0) = 0;

    % skip the thing if: see the
    if all(trial_tf(:)==0)
        disp('All values equal 0 after aperiodic subtraction');
        return;
    end
    
    % TF for iterating
    trial_tf_iter = trial_tf;

    while true
        % Compute noise floor
        thresh = 2 * std(trial_tf_iter(:));

        % Find peak
        [~,I] = max(trial_tf_iter(:));
        [peak_freq_idx, peak_time_idx] = ind2sub(size(trial_tf_iter),I);
        peak_freq = search_freqs(peak_freq_idx);
        peak_amp_iter = trial_tf_iter(peak_freq_idx, peak_time_idx);
        peak_amp_base = trial_tf(peak_freq_idx, peak_time_idx);
        if peak_amp_iter < thresh
            break
        end

        % Fit 2D Gaussian and subtract from TF
        [rloc, llec, uloc, dloc] = fwhm_burst_norm(trial_tf_iter,...
            [peak_freq_idx, peak_time_idx]);

        % REMOVE DEGENERATE GAUSSIAN
        vert_isnan = any(isnan([uloc, dloc]));
        horiz_isnan = any(isnan([rloc, llec]));
        if vert_isnan
            v_sh = round((length(search_freqs) - peak_freq_idx) / 2);
            if v_sh <= 0
                v_sh = 1;
            end
            uloc = v_sh;
            dloc = v_sh;

        elseif horiz_isnan
            h_sh = round((length(times) - peak_time_idx) / 2);
            if h_sh <= 0
                h_sh = 1;
            end
            rloc = h_sh;
            llec = h_sh;
        end

        hv_isnan = any([vert_isnan, horiz_isnan]);

        fwhm_f_idx = uloc + dloc;
        fwhm_f = (search_freqs(2)-search_freqs(1))*fwhm_f_idx;
        fwhm_t_idx = llec + rloc;
        fwhm_t = (times(2) - times(1))*fwhm_t_idx;
        sigma_t = (fwhm_t_idx) / 2.355;
        sigma_f = (fwhm_f_idx) / 2.355;
        z = peak_amp_iter * gaus2d(x_idx, y_idx, peak_time_idx,...
            peak_freq_idx, sigma_t, sigma_f);
        new_trial_TF_iter = trial_tf_iter - z;

        if peak_freq>=band_lims(1) && peak_freq<=band_lims(2) && ~hv_isnan
            % Extract raw burst signal
            dur = [max([1, peak_time_idx - llec]),...
                min([length(raw_trial), peak_time_idx + rloc])];
            raw_signal = raw_trial(dur(1):dur(2));

            % Bandpass filter
            freq_range = [max([1, peak_freq_idx - dloc]),...
                min([length(search_freqs) , peak_freq_idx + uloc])];
            
            dc=mean(raw_signal);
            % Pad with 1s on either side
            padded_data=[repmat(dc, 1, sfreq) raw_signal repmat(dc, 1, sfreq)];
            filtered = ft_preproc_bandpassfilter(padded_data, sfreq,...
                search_freqs(freq_range), 6, 'but', 'twopass', 'reduce');       
            filtered=filtered(sfreq+1:sfreq+length(raw_signal));
            
            % Hilbert transform
            analytic_signal = hilbert(filtered);
            % Get phase
            instantaneous_phase = mod(unwrap(angle(analytic_signal)), pi);
            
            % Find local phase minima with negative deflection closest to TF peak
            % If no minimum is found, the error is caught and no burst is added
            [~,zero_phase_pts]= findpeaks(-1*instantaneous_phase);
            if isempty(zero_phase_pts)
                adjustment=inf;
            else
                [~,min_idx]=min(abs((dur(2) - dur(1)) * .5 - zero_phase_pts));
                closest_pt = zero_phase_pts(min_idx);
                new_peak_time_idx = dur(1) + closest_pt;
                adjustment = (new_peak_time_idx - peak_time_idx) * 1 / sfreq;
            end

            % Keep if adjustment less than 30ms
            if abs(adjustment) < .03

                % If burst won't be cutoff
                if new_peak_time_idx > half_wlen && new_peak_time_idx + half_wlen < length(raw_trial)
                    peak_time = times(new_peak_time_idx);

                    overlapped=false;
                    % Check for overlap
                    for b_idx=1:length(bursts.peak_time)
                        o_t=bursts.peak_time(b_idx);
                        o_fwhm_t=bursts.fwhm_time(b_idx);
                        if overlap([peak_time-.5*fwhm_t, peak_time+.5*fwhm_t],...
                                [o_t-.5*o_fwhm_t, o_t+.5*o_fwhm_t])
                            overlapped=true;
                            break
                        end
                    end

                    if ~overlapped
                        % Get burst
                        burst = raw_trial(new_peak_time_idx - half_wlen:new_peak_time_idx + half_wlen);
                        % Remove DC offset
                        burst = burst - mean(burst);
                        bursts.waveform_times = times(new_peak_time_idx - half_wlen:new_peak_time_idx + half_wlen) - times(new_peak_time_idx);

                        % Flip if positive deflection
                        [~,peak_idxs]= findpeaks(filtered);                        
                        peak_dists = abs(peak_idxs - closest_pt);
                        [~,trough_idxs]= findpeaks(-1*filtered);                        
                        trough_dists = abs(trough_idxs - closest_pt);

                        polarity=0;
                        if isempty(trough_dists) || (~isempty(peak_dists) && min(peak_dists) < min(trough_dists))
                            burst = burst*-1.0;
                            polarity=1;
                        end

                        bursts.waveform(end+1,:)=burst;
                        bursts.peak_freq(end+1)=peak_freq;
                        bursts.peak_amp_iter(end+1)=peak_amp_iter;
                        bursts.peak_amp_base(end+1)=peak_amp_base;
                        bursts.peak_time(end+1)=peak_time;
                        bursts.peak_adjustment(end+1)=adjustment;
                        bursts.fwhm_freq(end+1)=fwhm_f;
                        bursts.fwhm_time(end+1)=fwhm_t;
                        bursts.polarity(end+1)=polarity;
                    end
                end
            end
        end

        trial_tf_iter = new_trial_TF_iter;
    end
end
