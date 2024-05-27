function subTc_filtered = post_processing_subject_timecourses(subTcs, Tr, band, cutoff_lim, display)
    % POST_PROCESSING_SUBJECT_TIMECOURSES applies detrending, filtering, and zscoring of fMRI time-course data.
    %
    % USAGE:
    %   subTc_filtered = post_processing_subject_timecourses(subTcs, Tr, band, cutoff_lim, display)
    %
    % INPUT:
    %   subTcs - A 3D array of time-course data (subjects x time points x components)
    %   Tr  - Repetition time (Tr) in seconds
    %   band   - A 1x2 vector specifying the bandpass frequency range [low high] in Hz
    %   cutoff_lim - A 1x2 vector specifying the passband cuttoffs [low high]
    %   display - (Optional) A boolean flag to plot post-processed time courses if true
    %
    % OUTPUT:
    %   subTc_filtered - A 3D array of post-processed fMRI time courses.

    % Input validation
    if nargin < 4
        error('Not enough input arguments. Provide subTcs, Tr, band, and cutoff_lim.');
    end
    if ~isnumeric(subTcs) || ndims(subTcs) ~= 3
        error('subTcs must be a 3D numeric array.');
    end
    if ~isnumeric(Tr) || ~isscalar(Tr) || Tr <= 0
        error('Repetition time must be a positive scalar.');
    end
    if ~isnumeric(band) || numel(band) ~= 2 || band(1) <= 0 || band(2) <= 0 || band(1) >= band(2)
        error('band must be a numeric vector with two positive elements [low high] such that low < high.');
    end
    if ~isnumeric(cutoff_lim) || numel(cutoff_lim) ~= 2 || cutoff_lim(1) <= 0 || cutoff_lim(2) <= 0 || ...
        cutoff_lim(1) >= cutoff_lim(2) || cutoff_lim(1) >= band(1) || cutoff_lim(2) <= band(2)
        error('cutoff_lim must be a numeric vector with two positive elements [low high] such that low < high.');
    end
    if nargin < 5
        display = false;
    end

    % Filter design parameters
    Fs = 1 / Tr;
    nyquist = Fs / 2;
    Wp = band / nyquist;
    Ws = [cutoff_lim(1), cutoff_lim(2)] / nyquist;
    Rp = 3;  % Passband ripple (dB)
    Rs = 30; % Stopband attenuation (dB)

    % Determine the order of the filter
    [n, Wn] = buttord(Wp, Ws, Rp, Rs);

    % Calculate filter coefficients
    [fkernB, fkernA] = butter(n, Wn);

    % Check filter stability
    if isstable(fkernB, fkernA)
        disp('Filter is stable');
    else
        error('Filter is NOT stable!');
    end

    % Initialize the output array
    subTc_filtered = zeros(size(subTcs));

    % Apply the filter to each subject and component
    for sub = 1:size(subTcs, 1)
        for comp = 1:size(subTcs, 3)
            % Extract the time-course data for the current subject and component
            tc = squeeze(subTcs(sub, :, comp));

            % Detrend, filter, and z-score the time-course data
            tc_detrended = detrend(tc);
            tc_filtered = filtfilt(fkernB, fkernA, tc_detrended);
            subTc_filtered(sub, :, comp) = zscore(tc_filtered);
        end
    end

    % Display post-processed time courses if the display flag is set to true
    if display
        figure;
        num_plots = min(5, size(subTcs, 1)); % Display up to 5 subjects
        for i = 1:num_plots
            subplot(num_plots, 1, i);
            plot(zscore(subTc_filtered(i, :, 1)), LineWidth=3); % Plot the first component for each subject
            title(sprintf('Subject %d, Component 1', i));
            xlabel('Time Points');
            ylabel('Z-scored Signal');
        end
    end
end
